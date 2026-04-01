# tsmi — Time-Series Manifold Integration

`tsmi.jl` is a Julia script that fits a **physics-based latent manifold** to a
climate or sea-level time series and then decomposes the signal via **Orthogonal
Matching Pursuit (OMP)** modulation.

---

## Overview

The model has two stages:

1. **Manifold construction** — builds a continuous latent signal *M(t)* by
   integrating a tidal forcing (sum of sinusoids at tidal periods) multiplied by
   an annual Gaussian comb (impulse train).  The comb supports a linear and
   quadratic temporal drift to capture slow geophysical changes.  A 40-year
   spin-up suppresses transients before the analysis window.

2. **OMP modulation** — represents the target time series *I(t)* as a sparse
   linear combination of complex exponentials `exp(i·κ·M(t))` (the *dictionary*),
   plus additive annual and semi-annual sinusoids.  Up to 18 atoms are selected
   by greedy matching pursuit with ridge regularisation.

Manifold parameters are optimised by a **Metropolis–Hastings random search**
with an adaptive temperature schedule and grid-search pre-conditioning for the
annual phase and initial condition.

---

## Files

| File | Description |
|---|---|
| `tsmi.jl` | Main Julia script |
| `ts.json` | Active parameter file (tidal amplitudes/phases, OMP coefficients, manifold settings) |
| `t.json` | Seed/template parameter file |
| `11final.dat` | Target time series — Warnemünde MSL monthly anomalies (decimal year, mm) |
| `11manifold.dat` | Empirical reference manifold used for visual comparison |
| `model_fit.csv` | Script output — columns: `t`, `I`, `I_model`, `Manifold`, `Manifold_Empirical` |
| `model_fit.png` | Three-panel diagnostic plot produced by `plot_model.py` |
| `plot_model.py` | Python visualisation script |

---

## Usage

```
julia tsmi.jl <timeseries.dat> <params.json> [options]
```

All options may be combined:

| Option | Effect |
|---|---|
| *(none)* | Full run: optimise manifold parameters, then run OMP on the target series |
| `--skip-optim` | Skip optimisation; load existing parameters and reconstruct the model |
| `--manifold` | Optimise the manifold only (no OMP step); saves the manifold time series to `m.dat` |
| `--final_mod <path>` | Skip manifold optimisation; load an external manifold file and run OMP directly |
| `--staged <cal.dat>` | Two-stage run: first fit the manifold to the calibration file `cal.dat` (stopping at *r* = 0.96), then apply the manifold to the main target series |

### Typical workflow

```bash
# 1. Optimise the manifold against the built-in target
julia tsmi.jl 11final.dat ts.json --manifold

# 2. Inspect m.dat, then run OMP modulation against the saved manifold
julia tsmi.jl 11final.dat ts.json --final_mod m.dat

# 3. Visualise the result
python plot_model.py        # writes model_fit.png
```

---

## Parameter file (`ts.json`)

The JSON file stores all model state and is read and updated in place on every
run (a timestamped `.bak_YYYYMMDD_HHMMSS` backup is written first).

Key fields:

| Field | Description |
|---|---|
| `tides` | Array of tidal components — each has `period` (days), `amplitude`, `phase` |
| `annual_phase` | Phase offset of the annual Gaussian impulse comb (years) |
| `primary_impulse_amp` | Amplitude of the primary annual impulse |
| `secondary_impulse_amp` | Amplitude of the half-year secondary impulse |
| `impulse_sigma` | Gaussian width of each impulse (years; default ≈ 1/78 ≈ one week) |
| `year_length` | Nominal year length used to convert tidal periods (days) |
| `drift_linear` | Linear drift of impulse timing relative to pivot time |
| `drift_quad` | Quadratic drift of impulse timing |
| `drift_pivot_time` | Reference time for the drift terms (set automatically on first run) |
| `ramp_slope` | Dissipative damping coefficient in the integrator |
| `dc_offset` | Constant offset added to the manifold |
| `initial_condition` | `{"A": …, "B": …}` — integrator state at `ic_ref_index` |
| `ic_ref_index` | Time index at which the initial condition is anchored |
| `initial_condition_ref_time` | Decimal-year time corresponding to `ic_ref_index` (used to re-anchor IC when the time series file changes) |
| `coefficients` | Sparse OMP coefficients — list of `{"k": κ, "re": …, "im": …}` |
| `additive_coefficients` | Annual/semi-annual additive term coefficients (`real`, `imag` arrays, length 4) |
| `comp_annual_amp` / `comp_annual_phase` | Compensatory annual sinusoid (added inside the integrator) |
| `comp_semi_amp` / `comp_semi_phase` | Compensatory semi-annual sinusoid |

---

## Output files

* **`model_fit.csv`** — written every run; contains the observed series, the
  OMP reconstruction, the parametric manifold, and the empirical manifold
  (if `11manifold1.dat` is present).
* **`m.dat`** — written in `--manifold` mode; tab-delimited `t  M(t)` suitable
  as input to `--final_mod`.
* **`ts.json`** — updated in place with the latest optimised parameters and
  coefficients.

---

## Dependencies

Julia packages (add once with `] add <Package>`):

```
DelimitedFiles  JSON3  Dates  LinearAlgebra  Statistics  Optim  DSP
```

Python packages required for `plot_model.py`:

```
matplotlib  numpy
```
