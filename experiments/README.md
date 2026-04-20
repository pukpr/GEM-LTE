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

## Annual forcing-feature model

`annual_forcing_feature_model.py` is the companion procedure that keeps the same
inputs and most of the same report artifacts, but replaces the monthly latent NN
with a fixed low-rank **annual increment** forcing-feature model. It is meant
for cases where year-to-year level shifts appear sharper than the smooth
residual NN can easily represent, and where `mean_forcing.dat` should act as the
primary annualized comparison feature.

This is now the main transition path for experiments that treat
`mean_forcing.dat` as the primary annualized forcing and use the year-scale
component as a structured correction for missing factors.

In addition to the shared reconstruction and latent artifacts, the annual
variant writes basis-determination reports for the `mean_forcing.dat` branch:

- `mean_forcing_basis_determination.csv`
- `mean_forcing_basis_determination_summary.csv`
- `mean_forcing_basis_coefficients.csv`
- `annual_site_loadings.csv`

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
python3 annual_forcing_feature_model.py \
  --config annual_forcing_feature_config.yml \
  --output-dir annual_forcing_feature_outputs
```

## Deterministic predictor

`lte_predict.py` is the companion predictor to `annual_forcing_feature_model.py`.
It uses the saved forcing-basis coefficients, annual correction parameters, and
`mean_forcing.dat` to generate a deterministic monthly time series for a target
index subdirectory. If the target site was part of the training fit, it uses the
exact saved parameters; otherwise it transfers them by geo-weighted proximity.
It saves a full-range comparison plot by default and can optionally display it.

Typical usage:

```bash
cd experiments
python3 lte_predict.py 11 \
  --model-dir annual_forcing_feature_outputs
```

Optional display:

```bash
python3 lte_predict.py 11 \
  --model-dir annual_forcing_feature_outputs \
  --show-plot
```

To compare against the actual observed series from column 3 of
`lte_results.csv` instead of the model series in column 2:

```bash
python3 lte_predict.py 234 \
  --model-dir annual_forcing_feature_outputs \
  --target-series actual
```
