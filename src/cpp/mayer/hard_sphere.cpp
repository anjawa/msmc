#include "hard_sphere.h"

HardSphere::HardSphere(double sigma) : sigma_(sigma) {}

double HardSphere::fij(double rij) const {
    return (rij <= sigma_) ? -1.0 : 0.0;
}
