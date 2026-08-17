#ifndef MSMC_METROPOLIS_HASTING_H
#define MSMC_METROPOLIS_HASTING_H

#include "../mayer/pair_potential.h"
#include <map>
#include <memory>
#include <vector>

class Model {
public:
    Model(std::shared_ptr<PairPotential> potential, std::map<int, double> virial);

    std::shared_ptr<PairPotential> potential;
    std::map<int, double> virial;
};

class MetropolisHasting {
public:
    MetropolisHasting(Model target_model, Model ref_model);

    double sample_virial(int virial_no, int num_samples, int warmup, double proposal_sigma, int seed);
    std::vector<std::vector<double>> sample_virial_overlap(int virial_no, int num_samples, int warmup,
                                                           int seed);
private:
    Model target_model_;
    Model ref_model_;
};

#endif
