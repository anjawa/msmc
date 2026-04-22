#include "pair_potential.h"

#include <cmath>
#include <stdexcept>

double PairPotential::prob_fij(double rij) const {
    return std::abs(fij(rij));
}

double PairPotential::compute_integrand(const double* pos, int virial_no) const {
    static constexpr double origin[3] = {0.0, 0.0, 0.0};
    const double* p1 = pos;
    const double* p2 = pos + 3;
    const double* p3 = pos + 6;

    if (virial_no == 2) {
        double r12 = dist(origin, p1);
        return fij(r12);
    } else if (virial_no == 3) {
        double r12 = dist(origin, p1);
        double r13 = dist(origin, p2);
        double r23 = dist(p1, p2);
        return fij(r12) * fij(r13) * fij(r23);
    } else if (virial_no == 4) {
        double r12 = dist(origin, p1);
        double r13 = dist(origin, p2);
        double r23 = dist(p1, p2);
        double r34 = dist(p2, p3);
        double r14 = dist(origin, p3);
        double r24 = dist(p1, p3);
        double gamma = fij(r12) * fij(r23) * fij(r34) * fij(r14);
        gamma *= (3.0 + 6.0 * fij(r13) +  fij(r24) * fij(r13));
        return gamma;
    }
    throw std::invalid_argument("Virial coefficient not implemented");
}

double PairPotential::dist(const double* a, const double* b) {
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
