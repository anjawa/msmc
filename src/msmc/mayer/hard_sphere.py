import msmc._core as _core


class HardSphere(_core.HardSphere):
    """Hard-sphere pair potential.

    Mayer-f function is -1 for r <= sigma, 0 otherwise.

    Parameters
    ----------
    sigma : float
        Hard-sphere diameter.
    """

    def __init__(self, sigma):
        super().__init__(sigma)
