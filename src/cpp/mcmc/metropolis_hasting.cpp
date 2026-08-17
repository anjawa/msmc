#include "metropolis_hasting.h"

#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace {

// Metropolis chain sampling pi = |gamma_own|
// gamma_other is kept up to date for the overlap function
class OverlapChain {
public:
    OverlapChain(const PairPotential* own, const PairPotential* other, int virial_no, int seed,
                 int tune_interval)
        : own_(own), other_(other), virial_no_(virial_no),
          n_dof_((virial_no - 1) * 3), rng_(seed), pos_(n_dof_), pos_trial_(n_dof_),
          tune_interval_(tune_interval) {
        // Chain begins inside the integrand's support; redraw if gamma == 0.
        for (int attempt = 0; attempt < 1000; attempt++) {
            for (int i = 0; i < n_dof_; i++) {
                pos_[i] = 0.3 * normal_(rng_);
            }
            gamma_own_ = own_->compute_integrand(pos_.data(), virial_no_);
            if (std::abs(gamma_own_) > 0.0) break;
        }
        if (std::abs(gamma_own_) == 0.0) {
            throw std::runtime_error("No start configuration with nonzero integrand found");
        }
        gamma_other_ = other_->compute_integrand(pos_.data(), virial_no_);
    }

    // tune=true only during warmup. Afterwards the step size is frozen
    void step(bool tune = false) {
        for (int i = 0; i < n_dof_; i++) {
            pos_trial_[i] = pos_[i] + step_sigma_ * normal_(rng_);
        }
        double gamma_new = own_->compute_integrand(pos_trial_.data(), virial_no_);
        double accept_prob = std::min(1.0, std::abs(gamma_new) / std::abs(gamma_own_));
        if (uniform_(rng_) < accept_prob) {
            pos_.swap(pos_trial_);
            gamma_own_ = gamma_new;
            gamma_other_ = other_->compute_integrand(pos_.data(), virial_no_);
            if (tune) {
                n_accepted_++;
            }
        }
        if (tune) {
            n_proposed_++;
        }
        if (tune && n_proposed_ == tune_interval_) {
            // Target 30-40 % acceptance.
            double rate = static_cast<double>(n_accepted_) / n_proposed_;
            if (rate < 0.30) {
                step_sigma_ /= 1.15;
            } else if (rate > 0.40) {
                step_sigma_ *= 1.15;
            }
            n_proposed_ = 0;
            n_accepted_ = 0;
        }
    }

    double gamma_own() const { return gamma_own_; }
    double gamma_other() const { return gamma_other_; }

private:
    const PairPotential* own_;
    const PairPotential* other_;
    int virial_no_;
    int n_dof_;
    std::mt19937 rng_;
    std::normal_distribution<double> normal_{0.0, 1.0};
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    std::vector<double> pos_;
    std::vector<double> pos_trial_;
    double gamma_own_ = 0.0;
    double gamma_other_ = 0.0;
    double step_sigma_ = 1.0;
    int n_proposed_ = 0;
    int n_accepted_ = 0;
    int tune_interval_;
};

// Sample means of one measurement run over both chains.
struct ChainAverages {
    double sign_t;
    double over_t;
    double sign_r;
    double over_r;
};

// Bennett overlap function g_OS = |g_t| |g_r| / (alpha |g_r| + |g_t|)
// (Benjamin, Schultz & Kofke 2007).
// Runs both chains for n_steps and returns the four sample means:
// <g_t/pi_T>, <g_OS/pi_T> (chain T) and <g_r/pi_R>, <g_OS/pi_R> (chain R).
ChainAverages run_chains(OverlapChain& chain_t, OverlapChain& chain_r, int n_steps, double alpha) {
    double sum_sign_t = 0.0, sum_over_t = 0.0;
    double sum_sign_r = 0.0, sum_over_r = 0.0;
    for (int n = 0; n < n_steps; n++) {
        chain_t.step();
        double abs_t = std::abs(chain_t.gamma_own());
        double abs_r = std::abs(chain_t.gamma_other());
        sum_sign_t += chain_t.gamma_own() / abs_t;
        sum_over_t += abs_r / (alpha * abs_r + abs_t);

        chain_r.step();
        abs_r = std::abs(chain_r.gamma_own());
        abs_t = std::abs(chain_r.gamma_other());
        sum_sign_r += chain_r.gamma_own() / abs_r;
        sum_over_r += abs_t / (alpha * abs_r + abs_t);
    }
    return {sum_sign_t / n_steps, sum_over_t / n_steps,
            sum_sign_r / n_steps, sum_over_r / n_steps};
}

}


Model::Model(std::shared_ptr<PairPotential> potential, std::map<int, double> virial)
    : potential(potential), virial(virial) {}

MetropolisHasting::MetropolisHasting(Model target_model, Model ref_model)
    : target_model_(target_model), ref_model_(ref_model) {}


double MetropolisHasting::sample_virial(int virial_no, int num_samples, int warmup, double proposal_sigma, int seed) {

    int n_particles = virial_no - 1;

    std::mt19937 rng(seed); // seed
    std::normal_distribution<double> normal(0.0, proposal_sigma);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::vector<double> positions(num_samples * n_particles * 3, 0.0);
    std::vector<double> gamma(num_samples, 0.0);
    std::vector<double> gamma_0(num_samples, 0.0);

    double* pos = positions.data();
    for (int i = 0; i < n_particles * 3; i++) {
        pos[i] = normal(rng);
    }
    gamma[0] = target_model_.potential->compute_integrand(pos, virial_no);
    gamma_0[0] = ref_model_.potential->compute_integrand(pos, virial_no);

    for (int n=1; n<num_samples; n++) {
        double* pos_old = positions.data() + (n-1) * n_particles * 3;
        double* pos_new = positions.data() + n * n_particles * 3;

        for (int i = 0; i < n_particles * 3; i++) {
            pos_new[i] = pos_old[i] + normal(rng);
        }

        double pdf_old = std::abs(target_model_.potential->compute_integrand(pos_old, virial_no));
        double pdf_new = std::abs(target_model_.potential->compute_integrand(pos_new, virial_no));

        double accept_prob = (pdf_old == 0.0) ? 1.0 : std::min(1.0, pdf_new / pdf_old);

        double uniform_randn = uniform(rng);

        if (uniform_randn < accept_prob) {
            gamma[n] = target_model_.potential->compute_integrand(pos_new, virial_no);
            gamma_0[n] = ref_model_.potential->compute_integrand(pos_new, virial_no);
        }
        else {
            std::copy(pos_old, pos_old + n_particles * 3, pos_new);
            gamma[n] = gamma[n - 1];
            gamma_0[n] = gamma_0[n - 1];
        }
    }

    double sum_ratio = 0.0;
    double sum_ratio_0 = 0.0;

    for (int i = warmup; i < num_samples; i++) {
        double pi = std::abs(gamma[i]);
        if (pi > 0.0) {
            sum_ratio += gamma[i] / pi;
            sum_ratio_0 += gamma_0[i] / pi;
        }
    }

    target_model_.virial[virial_no] = ref_model_.virial.at(virial_no) * sum_ratio / sum_ratio_0;
    return target_model_.virial.at(virial_no);
}

std::vector<std::vector<double>> MetropolisHasting::sample_virial_overlap(int virial_no, int num_samples,
                                                                          int warmup, int seed) {

    // Aim for approx. 100 tuning windows during warmup (min 50, max 200 proposals per window).
    int tune_interval = std::clamp(warmup / 100, 50, 200);

    OverlapChain chain_t(target_model_.potential.get(), ref_model_.potential.get(), virial_no, seed,
                         tune_interval);
    OverlapChain chain_r(ref_model_.potential.get(), target_model_.potential.get(), virial_no, seed + 1,
                         tune_interval);

    for (int n = 0; n < warmup; n++) {
        chain_t.step(true);
        chain_r.step(true);
    }

    // alpha from short pre-runs via Bennett's criterion
    // (Benjamin, Schultz & Kofke 2007): alpha = <g_OS/pi_R> / <g_OS/pi_T>.
    // alpha only needs to be roughly right, so the pre-runs are capped.
    double alpha = 1.0;
    int n_pre = std::max(1, std::min(num_samples / 20, 50000));
    for (int it = 0; it < 2; it++) {
        ChainAverages avg = run_chains(chain_t, chain_r, n_pre, alpha);
        double alpha_new = avg.over_r / avg.over_t;
        if (std::isfinite(alpha_new) && alpha_new > 0.0) {
            alpha = alpha_new;
        }
    }

    // batch means of the four averages
    int batch_len = std::max(1, num_samples / 10000);
    int n_batches = num_samples / batch_len;
    std::vector<double> batch_sign_t(n_batches), batch_over_t(n_batches);
    std::vector<double> batch_sign_r(n_batches), batch_over_r(n_batches);
    for (int b = 0; b < n_batches; b++) {
        ChainAverages avg = run_chains(chain_t, chain_r, batch_len, alpha);
        batch_sign_t[b] = avg.sign_t;
        batch_over_t[b] = avg.over_t;
        batch_sign_r[b] = avg.sign_r;
        batch_over_r[b] = avg.over_r;
    }

    return {batch_sign_t, batch_over_t, batch_sign_r, batch_over_r};
}
