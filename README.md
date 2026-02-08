# GEM-LTE

GeoEnergyMath Laplace's Tidal Equation modeling and fitting software.

## Essential operation (ENSO_OPT ➜ LT.EXE)

The main Ada driver is `src/enso_opt.adb`. It is built into `lt.exe` (or `lt` on
Linux/macOS) via the GNAT project file `lte.gpr`.

1. **Build**: `gprbuild -P lte.gpr enso_opt` (or `gprbuild -P lte.gpr`).
2. **Stage the executable**: copy `obj/enso_opt.exe` to `run/lt.exe`
   (see `update_exe.bat` or `update_exe.ps1`).
3. **Runtime flow**:
   - Reads environment configuration such as `NUMBER_OF_PROCESSORS`, `DLOD_DAT`,
     `TIMEOUT`, `DLOD_SCALE`, and `EXPECT`.
   - Loads the dLOD reference data (default `../dLOD3.dat`) and parameter files
     through the shared primitives (the `run/lt.exe.par*` and
     `run/lt.exe.resp*` files are examples).
   - Starts the solver tasks, periodically printing status and checking for
     console input.
   - Interactive console controls: `q`/`s` stop, `x` exits without saving, and
     `1`–`9` trigger stored solution outputs.
   - The `TIMEOUT` cycle stops the run when the elapsed time is exceeded.

For a more detailed walkthrough, see the [wiki page](wiki/Operation.md).

## GUI workflow (Windows command line)

The interactive GUI lives in `experiments/Feb2026/lte_gui.py` and is launched
from a Windows command prompt (Linux support is planned; currently the launcher
uses `cmd.exe`). It helps you select climate indices and mean sea level (MSL)
sites to model and simulate.

1. Run the GUI from the `experiments\Feb2026` level (the script expects to run
   at that level and searches subdirectories):
   `python experiments\Feb2026\lte_gui.py`.
2. Choose the root directory that contains index folders (for example,
   `experiments\Feb2026`). The GUI lists each directory except `locs`,
   `scripts`, and `rlr_data`.
3. Select an index directory. If the directory name matches an entry in
   `ID.yml`, the GUI shows the MSL site name/country metadata.
4. Configure the `TIMEOUT`, metric (CC/DTW), and training interval.
5. Click **Run lt** to launch `lt.exe` in a new console with environment
   variables (`METRIC`, `TIMEOUT`, `TRAIN_START`, `TRAIN_STOP`, `CLIMATE_INDEX`,
   `IDATE`) set for the run.
6. Click **Run plot** to execute `plot.py` for the selected index and
   **Refresh PNG** to preview the newest plot in the GUI.

> Screenshot placeholder: insert GUI screenshot here.

## Release notes / packaging

Releases should package the compiled `lt.exe`, the `run/` data/parameter files,
and the experiment directories used with `lte_gui.py`. If you want a new GitHub
Release, tag the commit and create a release from the tag in GitHub.
