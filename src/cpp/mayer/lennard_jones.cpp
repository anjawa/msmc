#include "lennard_jones.h"

#include <cmath>

LennardJones::LennardJones(double sigma, double epsilon, double k, double T)
    : sigma_(sigma), epsilon_(epsilon), k_(k), T_(T) {}

double LennardJones::uij(double rij) const {
    double ratio = sigma_ / rij;
    double ratio6 = std::pow(ratio,6);
    return 4.0 * epsilon_ * (ratio6 * ratio6 - ratio6);
}

double LennardJones::fij(double rij) const {
    return std::exp(-uij(rij) / (k_ * T_)) - 1.0;
}
