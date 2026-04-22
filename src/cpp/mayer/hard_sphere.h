#ifndef MSMC_HARD_SPHERE_H
#define MSMC_HARD_SPHERE_H

#include "pair_potential.h"

class HardSphere : public PairPotential {
public:
    explicit HardSphere(double sigma);

    double fij(double rij) const override;

private:
    double sigma_;
};

#endif
