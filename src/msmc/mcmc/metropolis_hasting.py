import numpy as np

import msmc._core as _core


def standard_error(x):
    """
    Standard error of the mean of a correlated series.

    The autocovariances are computed via FFT; the integrated
    autocorrelation time tau is estimated with Geyer's initial monotone
    sequence (as in Stan's ESS diagnostic) and floored at 1.
    """
    x = np.asarray(x, dtype=float)
    m = x.size
    dev = x - x.mean()

    # autocovariances c(0..m-1); rfft zero-pads to 2m
    spectrum = np.fft.rfft(dev, 2 * m)
    c = np.fft.irfft(np.abs(spectrum) ** 2, 2 * m)[:m] / m

    if c[0] == 0.0:
        return 0.0  # constant series

    # tau = 1 + 2 sum_k rho(k), summed over pair sums rho(2k) + rho(2k+1)
    # while they stay positive and monotonically decreasing.
    n_pairs = m // 2
    pairs = (c[0 : 2 * n_pairs : 2] + c[1 : 2 * n_pairs : 2]) / c[0]
    negative = np.nonzero(pairs < 0)[0]
    if negative.size:
        pairs = pairs[: negative[0]]
    pairs = np.minimum.accumulate(pairs)
    tau = max(-1.0 + 2.0 * pairs.sum(), 1.0)

    return np.sqrt(c[0] * tau / m)


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
        self._ref_model = ref_model

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

    def sample_virial_overlap(
        self, virial_no=2, num_samples=100000, warmup=20000, seed=42
    ):
        """
        Overlap sampling MSMC (Benjamin, Schultz & Kofke, 2007).

        The proposal step size is tuned during warmup (target 30-40 %
        acceptance) and frozen afterwards.

        The standard error accounts for autocorrelation (as in Stan's ESS
        diagnostic).

        Parameters
        ----------
        virial_no : int, optional
            Virial coefficient number (default 2).
        num_samples : int, optional
            Number of production MCMC samples per chain (default 100000).
        warmup : int, optional
            Burn-in samples used for equilibration and step-size tuning
            (default 20000).
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        tuple[float, float]
            Estimated virial coefficient and its standard error.
        """
        sign_t, over_t, sign_r, over_r = super().sample_virial_overlap(
            virial_no, num_samples, warmup, seed
        )
        sign_t = np.asarray(sign_t)
        over_t = np.asarray(over_t)
        sign_r = np.asarray(sign_r)
        over_r = np.asarray(over_r)

        ratio = (sign_t.mean() / over_t.mean()) * (over_r.mean() / sign_r.mean())
        value = self._ref_model.virial[virial_no] * ratio

        # Error via linearization of the ratio f = s_t o_r / (o_t s_r)
        # we use: df/f = ds_t/s_t - do_t/o_t + do_r/o_r - ds_r/s_r
        influence = (
            sign_t / sign_t.mean()
            - over_t / over_t.mean()
            + over_r / over_r.mean()
            - sign_r / sign_r.mean()
        )
        error = abs(value) * standard_error(influence)

        return value, error

    @staticmethod
    def compute_mc_stats(x):
        """Compute mean and standard error of the mean."""
        mean = np.mean(x)
        sem = np.std(x, ddof=1) / np.sqrt(len(x))
        return [mean, sem]
