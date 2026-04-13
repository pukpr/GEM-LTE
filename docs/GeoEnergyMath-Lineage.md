# GeoEnergyMath Lineage

GEM-LTE evolved directly from the earlier
[GeoEnergyMath](https://github.com/pukpr/GeoEnergyMath) repository, which is now
deprecated and points here as its successor.  Both projects implement Laplace's
Tidal Equation (LTE) modelling in Ada/GNAT and share a common design lineage
described below.

---

## What carried over

### Core Ada source files

These files originated in GeoEnergyMath and continue in GEM-LTE (often with
enhancements):

- [`src/gem.ads` / `src/gem.adb`](../src/gem.ads) – root GEM package;
  utility types, constants, and helper subprograms
- [`src/gem-lte.ads` / `src/gem-lte.adb`](../src/gem-lte.ads) – top-level LTE
  package; tidal-period enumeration and constituent arrays
- [`src/gem-lte-primitives.ads` / `.adb`](../src/gem-lte-primitives.ads) –
  core algorithms: `Tide_Sum`, `LTE`, `IIR`/`FIR` filters, and every
  fit-quality metric (CC, DTW, EMD, …)
- [`src/gem-lte-primitives-shared.ads` / `.adb`](../src/gem-lte-primitives-shared.ads) –
  shared global state: parameter tables, response configuration, and the
  file I/O that reads `.par` / `.resp` files
- [`src/gem-lte-primitives-solution.ads` / `.adb`](../src/gem-lte-primitives-solution.ads) –
  solver loop that drives parameter optimisation across concurrent tasks
- [`src/gem-dlod.ads` / `src/gem-dlod.adb`](../src/gem-dlod.ads) – loads the
  dLOD (delta Length-of-Day) reference data that calibrates the tidal forcing
- [`src/gem-random_descent.ads` / `.adb`](../src/gem-random_descent.ads) –
  random-descent optimiser used to search the parameter space
- [`src/enso_opt.adb`](../src/enso_opt.adb) – main Ada driver (entry point);
  reads environment variables, starts solver tasks, and handles interactive
  console controls (`q`, `x`, `1`–`9`)
- [`lte.gpr`](../lte.gpr) – GNAT project file for building the `enso_opt`
  executable

### Reference data

- [`dLOD3.dat`](../dLOD3.dat) – the canonical dLOD time-series file;
  its earlier version (`dLOD.dat`) shipped with GeoEnergyMath and the same
  dataset underpins all models in GEM-LTE

### Design patterns

- **`.par` / `.resp` file convention** – parameter and response configuration
  files read at runtime, pioneered in GeoEnergyMath and preserved in
  [`run/`](../run/)
- **Interactive console controls** – keyboard shortcuts (`q`/`s`, `x`,
  `1`–`9`) described in both READMEs
- **Climate-index directory layout** – each index (ENSO, QBO, NAO, SOI, …)
  gets its own subdirectory containing `lt.exe`, `.par`, `.resp`, and plot
  output; inherited from GeoEnergyMath's per-index subdirectories and
  formalised in [`experiments/Feb2026/`](../experiments/Feb2026/)

---

## What GEM-LTE added or replaced

| GeoEnergyMath | GEM-LTE equivalent / improvement |
|---|---|
| `gem-lte-lib.ads` – hard-coded LTE constant arrays (QBO, LOD, gravity, …) | Replaced by runtime `.par` files; constants removed from source |
| `gem-ephemeris.adb/.ads` – planetary ephemeris helper | Not carried over; tidal periods are now fully parameterised in `.par` files |
| `gem-matrices.adb/.ads` – matrix utilities | Not carried over; regression uses the new `gem-mix_regression` package |
| Multiple single-purpose drivers (`all_opt`, `lod_opt`, `nao_opt`, `qbo_opt`, `tidal_opt`) | Unified into the single `enso_opt` driver with runtime index selection |
| Manual copy of the executable into each index directory | [`scripts/update_exe.ps1`](../scripts/update_exe.ps1) / [`setup_with_gnat.ps1`](../setup_with_gnat.ps1) automate staging |
| No GUI | Python GUI [`experiments/Feb2026/lte_gui.py`](../experiments/Feb2026/) for index/site selection and plot preview |
| `src/gem-lte-primitives.adb` – regression baked in | [`src/gem-mix_regression.adb/.ads`](../src/gem-mix_regression.ads) extracted as a reusable package |
| – | [`src/gem-lte-primitives-param_b_overlay.adb/.ads`](../src/gem-lte-primitives-param_b_overlay.ads) – overlay mechanism for parameter-B variants |
| – | [`src/index_regress.adb`](../src/index_regress.adb) – standalone index regression driver |
| – | [`TidalFactorBarChart.py`](../TidalFactorBarChart.py) – visualises tidal-factor amplitudes as a bar chart |

---

## See also

- [GeoEnergyMath repository](https://github.com/pukpr/GeoEnergyMath) (archived / deprecated)
- [Operation wiki page](Operation.md) – detailed walkthrough of the build and runtime workflow
- [GEM-LTE README](../README.md)
