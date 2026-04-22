import msmc._core as _core


class LennardJones(_core.LennardJones):
    """Lennard-Jones pair potential.

    Parameters
    ----------
    sigma : float
        Length scale (distance at which potential is zero).
    epsilon : float
        Energy scale (depth of potential well).
    k : float
        Boltzmann constant.
    T : float
        Temperature.
    """

    def __init__(self, sigma, epsilon, k, T):
        super().__init__(sigma, epsilon, k, T)

    def uij(self, rij):
        """
        Lennard-Jones pair energy.

        Parameters
        ----------
        rij : float
            Separation distance between molecule i and j.

        Returns
        -------
        float
            Pair potential energy.
        """
        return super().uij(rij)
