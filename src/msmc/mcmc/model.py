import msmc._core as _core


class Model(_core.Model):
    """
    Bundles a PairPotential with its virial coefficients.

    Parameters
    ----------
    potential : PairPotential
        Pair potential model (e.g. HardSphere, LennardJones).
    virial : dict[int, float], optional
        Known virial coefficients keyed by virial number (e.g. {2: b2, 3: b3}).
    """

    pass
