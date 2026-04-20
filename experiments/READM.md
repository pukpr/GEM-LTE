# Experiments

This directory contains exploratory numerical experiments that sit alongside the
main GEM-LTE workflows.

## Latent-manifold experiment

`nn_latent_manifold_experiment.py` trains a lightweight geo-aware neural model
on `Feb2026/*/lte_results.csv` time series to learn a shared latent manifold
across local sites and regional composites.

Current implementation highlights:

- uses `Feb2026/ID.yml` for site metadata and basin inference
- excludes configurable site prefixes and specific problematic sites
- supports `tanh`, `sin`, and `hybrid` temporal encoders
- supports `mean_forcing.dat` as either a hint or a stronger
  `composite_residual` baseline branch
- writes reconstruction metrics, latent trajectories, spectral diagnostics, and
  latent-vs-mean-forcing comparison outputs
- reports the fractional mean-correlation breakdown between the
  `mean_forcing.dat` branch and the residual NN branch for composite runs

Typical usage from the repository root:

```bash
cd experiments
python3 nn_latent_manifold_experiment.py \
  --config nn_latent_manifold_config.yml \
  --output-dir nn_latent_manifold_outputs
```

The default configuration lives in `nn_latent_manifold_config.yml`, and sweep
support is available through `nn_latent_manifold_sweep.py` with
`nn_latent_manifold_sweep.yml`.

## Annualized offshoot

`nn_latent_manifold_annual_experiment.py` is a companion procedure that keeps
the same inputs and most of the same report artifacts, but replaces the monthly
latent NN with a low-rank **annual increment** latent model. It is meant for
cases where year-to-year level shifts appear sharper than the smooth residual NN
can easily represent.

This is now the main transition path for experiments that treat
`mean_forcing.dat` as the primary annualized forcing and use the year-scale
component as a structured correction for missing factors.

In addition to the shared reconstruction and latent artifacts, the annual
variant writes basis-determination reports for the `mean_forcing.dat` branch:

- `mean_forcing_basis_determination.csv`
- `mean_forcing_basis_determination_summary.csv`

These report per-basis single-term `R^2` and drop-one `delta R^2` for the basis
library `{1, M, M^2, M^3, sin(M), cos(M), M sin(M), M cos(M)}` plus any
configured higher-harmonic forcing terms such as `sin(k M)` / `cos(k M)`,
alongside the separate annual residual-component metrics.

Higher harmonics can be supplied explicitly with
`mean_forcing_hint.extra_harmonic_factors`, or generated automatically with
`mean_forcing_hint.max_harmonic_factor`, which adds every harmonic from `k=2`
through that maximum. To discourage overfitting in the expanded forcing basis,
`mean_forcing_hint.lasso_alpha` applies an L1 sparsity penalty on all
non-constant basis terms.

Typical usage:

```bash
cd experiments
python3 nn_latent_manifold_annual_experiment.py \
  --config nn_latent_manifold_annual_config.yml \
  --output-dir nn_latent_manifold_annual_outputs
```
