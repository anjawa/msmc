import argparse

import numpy as np

from msmc.mayer import LennardJones, HardSphere
from msmc.mcmc import MetropolisHasting, Model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--T_reduced", type=float, default=1.0)
    args = parser.parse_args()

    sigma_lj = 1.0
    sigma_hs = sigma_lj

    # Reduced units: epsilon=1, k=1, T=T_reduced
    lj = LennardJones(sigma=sigma_lj, epsilon=1.0, k=1.0, T=args.T_reduced)
    hs = HardSphere(sigma=sigma_hs)

    # HS reference virial coefficients (Ree & Hoover, 1966)
    b = 2.0 * np.pi / 3.0 * sigma_hs**3
    b_ref = {
        2: b,
        3: 5.0 / 8.0 * b**2,
        4: (
            -89.0 / 280.0
            + (219.0 * np.sqrt(2) + 4131.0 * np.arccos(1.0 / np.sqrt(3)))
            / 2240.0
            / np.pi
        )
        * b**3,
        5: 0.1103 * b**4,
        6: 0.0386 * b**5,
    }

    target = Model(potential=lj, virial={})
    ref = Model(potential=hs, virial=b_ref)

    msmc = MetropolisHasting(target_model=target, ref_model=ref)

    b2 = msmc.sample_virial(
        virial_no=2, num_samples=int(1e6), proposal_sigma=1.0, seed=42
    )
    b3 = msmc.sample_virial(
        virial_no=3, num_samples=int(1e6), proposal_sigma=1.0, seed=42
    )
    b4 = msmc.sample_virial(
        virial_no=4, num_samples=int(1e6), proposal_sigma=1.0, seed=42
    )

    print("B2 = {}; B2/b = {}".format(b2, b2 / b))
    print("B3 = {}; B3/b^2 = {}".format(b3, b3 / b**2))
    print("B4 = {}; B4/b^3 = {}".format(b4, b4 / b**3))
