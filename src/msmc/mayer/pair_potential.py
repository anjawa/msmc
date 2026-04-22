import numpy as np

import msmc._core as _core


class PairPotential(_core.PairPotential):
    """Base class for pairwise-additive potential models."""

    def fij(self, rij):
        """
        Mayer-f function for a pair of molecules.

        Parameters
        ----------
        rij : float
            Separation distance between molecule i and j.

        Returns
        -------
        float
            Mayer-f function value.
        """
        return super().fij(rij)

    def prob_fij(self, rij):
        """
        Sampling probability for a pair of molecules.

        Parameters
        ----------
        rij : float
            Separation distance between molecule i and j.

        Returns
        -------
        float
            Absolute value of Mayer-f function.
        """
        return super().prob_fij(rij)

    def compute_integrand(self, pos, virial_no):
        """
        Compute cluster integral integrand for virial coefficients B2-B4.

        Cluster integral corresponds to Mayer-f function for
        pairwise-additive models (Benjamin et al, 2007).

        Parameters
        ----------
        pos : np.ndarray, shape (N, 3)
            Position vectors of particles relative to particle 1 at origin.
        virial_no : int
            Virial coefficient number. 2 <= virial_no <= 4.

        Returns
        -------
        float
            Integrand value of cluster integral.
        """
        return super().compute_integrand(
            np.ascontiguousarray(pos, dtype=np.float64), virial_no
        )
