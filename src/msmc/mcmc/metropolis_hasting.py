import numpy as np

import msmc._core as _core


class MetropolisHasting(_core.MetropolisHasting):
    """
    Mayer Sampling Monte Carlo (Singh & Kofke, 2004).

    Computes virial coefficients using Metropolis-Hastings MCMC
    with a hard-sphere reference system.

    Parameters
    ----------
    target_model : Model
        Target model (PairPotential + virial coefficients).
    ref_model : Model
        Reference model (PairPotential + virial coefficients).
    """

    def __init__(self, target_model, ref_model):
        super().__init__(target_model=target_model, ref_model=ref_model)

    def sample_virial(
        self, virial_no=2, num_samples=10000, warmup=200, proposal_sigma=0.7, seed=42
    ):
        """
        Direct sampling MSMC (Singh & Kofke, 2004).

        Parameters
        ----------
        virial_no : int, optional
            Virial coefficient number (default 2).
        num_samples : int, optional
            Number of MCMC samples (default 10000).
        warmup : int, optional
            Burn-in samples to discard (default 200).
        proposal_sigma : float, optional
            Std of Gaussian proposal distribution. Tune so that
            acceptance rate is ~50% (Singh & Kofke, 2004).
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        float
            Estimated virial coefficient.
        """
        return super().sample_virial(
            virial_no, num_samples, warmup, proposal_sigma, seed
        )

    @staticmethod
    def compute_mc_stats(x):
        """Compute mean and standard error of the mean."""
        mean = np.mean(x)
        sem = np.std(x, ddof=1) / np.sqrt(len(x))
        return [mean, sem]
