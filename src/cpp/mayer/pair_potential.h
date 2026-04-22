#ifndef MSMC_PAIR_POTENTIAL_H
#define MSMC_PAIR_POTENTIAL_H

class PairPotential {
public:
    virtual ~PairPotential() = default;

    virtual double fij(double rij) const = 0;
    virtual double prob_fij(double rij) const;

    double compute_integrand(const double* pos, int virial_no) const;

private:
    static double dist(const double* a, const double* b);
};

#endif
