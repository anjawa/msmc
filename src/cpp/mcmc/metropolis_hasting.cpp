#include "metropolis_hasting.h"

#include <cmath>
#include <random>
#include <vector>
#include <algorithm>


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
