#ifndef MSMC_LENNARD_JONES_H
#define MSMC_LENNARD_JONES_H

#include "pair_potential.h"

class LennardJones : public PairPotential {
public:
    LennardJones(double sigma, double epsilon, double k, double T);

    double uij(double rij) const;
    double fij(double rij) const override;

private:
    double sigma_;
    double epsilon_;
    double k_;
    double T_;
};

#endif
