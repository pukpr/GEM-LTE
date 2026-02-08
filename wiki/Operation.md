# GEM-LTE Operation Overview

This page documents the essential flow of the GEM-LTE software, starting from
the Ada driver (`ENSO_OPT.adb`) and ending with the interactive GUI launched
from a Windows command line.

## 1) Build chain: ENSO_OPT.adb ➜ LT.EXE

1. **Compile** with GNAT:
   - `gprbuild -P lte.gpr enso_opt`
2. **Stage the binary**:
   - Copy `obj/enso_opt.exe` to `run/lt.exe` (see `update_exe.bat` or
     `update_exe.ps1`).

The `lte.gpr` project lists `enso_opt.adb` as one of the main programs, so the
resulting executable is `enso_opt.exe`, which is renamed to `lt.exe` for use by
the GUI and scripts.

## 2) Runtime behavior of LT.EXE

When `lt.exe` starts (from `ENSO_OPT.adb`), it:

- Reads environment variables for configuration.
- Loads the dLOD reference data (`DLOD_DAT`, default `../dLOD3.dat`).
- Loads shared parameter/response state if present (see `run/lt.exe.par*` and
  `run/lt.exe.resp*` for examples).
- Starts solver tasks and periodically prints status updates.
- Listens for interactive console input.

### Key environment variables

| Variable | Purpose |
| --- | --- |
| `NUMBER_OF_PROCESSORS` | Controls solver parallelism. |
| `DLOD_DAT` | Path to the dLOD reference file (default `../dLOD3.dat`). |
| `TIMEOUT` | Maximum runtime before the solver is stopped. |
| `DLOD_SCALE` | Scale factor applied to dLOD amplitudes. |
| `EXPECT` | Uses blocking input instead of immediate input when `true`. |

### Interactive console controls

| Key | Action |
| --- | --- |
| `q` / `s` | Stop the solver. |
| `x` | Exit immediately (no save). |
| `1`–`9` | Trigger stored solution outputs. |

## 3) GUI workflow (Windows command line)

The GUI in `experiments/Feb2026/lte_gui.py` launches `lt.exe` and helps select
climate indices and MSL sites.

1. Launch from **Windows Command Prompt** at the `experiments\Feb2026` level
   (the script expects to run from there and searches subdirectories):
   - `python .\lte_gui.py`
2. Choose a root folder containing index directories (e.g.,
   `experiments\Feb2026`). The GUI lists all child directories except
   `locs`, `scripts`, and `rlr_data`.
3. Select an index directory. If its name matches an entry in `ID.yml`, the GUI
   shows the MSL site name and country.
4. Configure:
   - `TIMEOUT` (slider)
   - Metric (`CC` or `DTW`)
   - Training interval start/stop
5. Click **Run lt** to spawn `lt.exe` in a new console window with environment
   variables set (`METRIC`, `TIMEOUT`, `TRAIN_START`, `TRAIN_STOP`,
   `CLIMATE_INDEX`, `IDATE`).
6. Click **Run plot** to execute `plot.py` for the selected index and
   **Refresh PNG** to preview the latest output in the GUI.

> Screenshot placeholder: add a GUI screenshot here.

### Linux compatibility note

The GUI currently spawns `cmd.exe`. To run on Linux, replace the command
launcher with a shell/terminal alternative that can open a new console window.

## 4) Release packaging guidance

When creating a new release, include:

- The compiled `lt.exe` (or `lt`).
- The `run/` parameter/reference files.
- The experiment folders (for example `experiments/Feb2026`).

Create a GitHub Release from a tagged commit so users can download a consistent
set of binaries and data.
