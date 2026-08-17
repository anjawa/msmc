"""Validate LJ virial coefficients (T* = 1) against Barker, Leonard & Pompe (1966)."""

import numpy as np
import pytest

from msmc.mayer import HardSphere, LennardJones
from msmc.mcmc import MetropolisHasting, Model

SIGMA = 1.0
B = 2.0 * np.pi / 3.0 * SIGMA**3

# Hard-sphere reference virial coefficients
B_REF = {
    2: B,
    3: 5.0 / 8.0 * B**2,
    4: 0.28695 * B**3,
}

# Barker, Leonard & Pompe, J. Chem. Phys. 44, 4206 (1966), T* = 1
LITERATURE = {
    2: -2.5381 * B,
    3: 0.4297 * B**2,
    4: -0.2769 * B**3,
}


@pytest.fixture(scope="module")
def msmc():
    lj = LennardJones(sigma=SIGMA, epsilon=1.0, k=1.0, T=1.0)
    hs = HardSphere(sigma=SIGMA)
    target = Model(potential=lj, virial={})
    ref = Model(potential=hs, virial=B_REF)
    return MetropolisHasting(target_model=target, ref_model=ref)


@pytest.mark.parametrize(
    "virial_no, num_samples, seeds",
    [
        (2, 1_000_000, (1, 2, 3)),
        (3, 1_000_000, (1, 2, 3)),
        (4, 5_000_000, (1, 2, 3)),
    ],
)
def test_overlap_sampling_vs_literature(msmc, virial_no, num_samples, seeds):
    results = [
        msmc.sample_virial_overlap(
            virial_no=virial_no, num_samples=num_samples, warmup=20000, seed=seed
        )
        for seed in seeds
    ]
    values = np.array([v for v, _ in results])
    errors = np.array([e for _, e in results])

    value = values.mean()
    error = np.sqrt(np.sum(errors**2)) / len(seeds)

    # Error must be small enough for the 3-sigma check to be meaningful.
    assert error < 0.35 * abs(LITERATURE[virial_no])
    assert abs(value - LITERATURE[virial_no]) < 3.0 * error
