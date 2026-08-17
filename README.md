# Mayer-sampling Monte Carlo

Mayer-sampling Monte Carlo for estimating virial coefficients.

This project is in early development.

## Background

This implementation is based on the Mayer sampling method introduced by Singh and Kofke:

> J. K. Singh and D. A. Kofke, "Mayer Sampling: Calculation of Cluster Integrals using Free-Energy Perturbation Methods," Phys. Rev. Lett. 92, 220601 (2004).

The recommended estimator uses overlap sampling, for details see:

> K. M. Benjamin, A. J. Schultz, and D. A. Kofke, "Higher-Order Virial Coefficients of Water Models," J. Phys. Chem. B 111, 11463 (2007).

Hard-sphere reference virial coefficients are taken from:

> F. H. Ree and W. G. Hoover, "Seventh Virial Coefficients for Hard Spheres and Hard Disks," J. Chem. Phys. 46, 4181–4197 (1967).

Compared to the original method, this implementation:
- Exploits symmetry in the cluster integrands instead of full explicit summation over all graphs.
- Perturbs all particles simultaneously per MC step, rather than a random subset.
- Uses the full cluster sum (rather than just the ring cluster) as the hard-sphere reference integrand, with known HS virial coefficients from Ree & Hoover as normalization.

Two estimators are available:
- `sample_virial_overlap` (recommended): overlap sampling (Benjamin, Schultz & Kofke, 2007). Step size is tuned during warmup.
- `sample_virial`: direct sampling of the target integrand (Singh & Kofke), fixed step size.

## Results vs Literature
The results were compared with data from:
> J. A. Barker, P. J. Leonard, and A. Pompe, J. Chem. Phys. 44, 4206 (1966)
> (as cited in Singh and Kofke)

Overlap sampling, $T^* = 1$, seed 42 ($10^6$ samples for $B_2$/$B_3$, $10^7$ for $B_4$; output of `scripts/msmc_example.py`):

| Virial coefficient | Literature | Implementation |
|:---|---:|---:|
| $B_2 / b$ | -2.5381 | -2.516 ± 0.023 |
| $B_3 / b^2$ | 0.4297 | 0.421 ± 0.024 |
| $B_4 / b^3$ | -0.2769 | -0.207 ± 0.092 |

## Getting started

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```bash
uv sync
```

Run the example script:

```bash
uv run python scripts/msmc_example.py --T_reduced 1.0
```

## Development

Set up pre-commit hooks:

```bash
uv run pre-commit install
```

For clangd support, install [just](https://github.com/casey/just) and run `just configure`.
