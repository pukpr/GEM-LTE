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
