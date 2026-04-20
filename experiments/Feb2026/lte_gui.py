#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import sys
import json
import shlex
import shutil
import subprocess
import threading
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime


# ---------------------------------------------------------------------------
# Core long period tides computation
# ---------------------------------------------------------------------------

def load_lpap(target_dir: str) -> list:
    """
    Load the 'lpap' list from 'lt.exe.p' in *target_dir*.

    The file is expected to be valid JSON.  Two formats are accepted:
      • A JSON object with a top-level key "lpap" whose value is the list.
      • A bare JSON array (the list itself).

    Each element of the list must be indexable as:
        element[0] – period in days
        element[1] – amplitude
        element[2] – phase in radians

    Returns
    -------
    list of [period_days, amplitude, phase] sublists / tuples
    """
    path = os.path.join(target_dir, "lt.exe.p")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"'lt.exe.p' not found in:\n  {target_dir}")

    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    if isinstance(data, dict):
        if "lpap" not in data:
            raise KeyError("Key 'lpap' not found in lt.exe.p JSON object.")
        lpap = data["lpap"]
    elif isinstance(data, list):
        lpap = data
    else:
        raise TypeError("lt.exe.p must contain a JSON object or a JSON array.")

    if not lpap:
        raise ValueError("The 'lpap' list is empty.")

    return lpap


def compose_sinusoids(lpap: list, t_days: np.ndarray) -> np.ndarray:
    """
    y(t) = Σ_i  amplitude_i · sin( 2π / period_i · t  +  phase_i )

    Parameters
    ----------
    lpap    : list of [period_days, amplitude, phase_radians]
    t_days  : 1-D array of time values in days

    Returns
    -------
    y : ndarray of the same shape as *t_days*
    """
    y = np.zeros_like(t_days, dtype=float)
    for entry in lpap:
        period_days = float(entry[0])
        amplitude   = float(entry[1])
        phase       = float(entry[2])
        if period_days == 0.0:
            continue
        omega = 2.0 * np.pi / period_days
        y += amplitude * np.sin(omega * t_days + phase)
    return y


# ---------------------------------------------------------------------------
# Plot helper
# ---------------------------------------------------------------------------

DAYS_PER_YEAR = 365.25  # mean Julian year (accounts for leap years)


def _date_array(start_year: int, n_years: int, dt_days: float = 1.0):
    """Return (t_days, dates) arrays for *n_years* at *dt_days* resolution."""
    t0 = datetime.date(start_year, 1, 1)
    n_steps = int(n_years * DAYS_PER_YEAR / dt_days)
    t_days = np.arange(n_steps, dtype=float) * dt_days
    dates = [t0 + datetime.timedelta(days=float(d)) for d in t_days]
    return t_days, dates


def plot_lpap(lpap: list, index: str, start_year: int = 1950, n_years: int = 50):
    """
    Plot the composed sum of sinusoids defined by *lpap* over *n_years*
    starting at *start_year*.
    """
    t_days, dates = _date_array(start_year, n_years, dt_days=1.0)
    y = compose_sinusoids(lpap, t_days)

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(dates, y, lw=0.8, color="steelblue")
    ax.axhline(0, color="k", lw=0.4)
    idx = Path(index).stem
    ax.set_title(
        f"{idx} — composed sum of {len(lpap)} tidal factors\n"
        f"{start_year} – {start_year + n_years}",
        fontsize=11,
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Amplitude (a.u.)")
    ax.xaxis.set_major_locator(mdates.YearLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    plt.xticks(rotation=45, ha="right")
    ax.grid(True, ls=":", lw=0.4, alpha=0.6)
    fig.tight_layout()
    plt.show()


def plot_lpap_amplitudes(lpap: list, index: str):
    """
    Horizontal bar chart of absolute-value amplitudes for each tidal
    periodicity in *lpap*.  Phase is excluded.

    Each bar is labelled with the period (in days) on the y-axis and the
    x-axis tick labels are formatted to 5 significant digits.
    """
    periods = [float(e[0]) for e in lpap]
    amplitudes = [abs(float(e[1])) for e in lpap]

    # Build y-axis labels: "period" in days, formatted to 5 sig-figs
    labels = [f"{p:.5g} d" for p in periods]

    n = len(amplitudes)
    y_pos = np.arange(n)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * n)))
    ax.barh(y_pos, amplitudes, align="center", color="steelblue", edgecolor="none")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()  # longest period at the top

    # Format x-axis ticks to 5 significant digits
    def _sig5(x, _pos):
        if x == 0:
            return "0"
        from math import log10, floor
        digits = 5 - 1 - floor(log10(abs(x)))
        digits = max(0, digits)
        return f"{x:.{digits}f}"

    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(_sig5))

    ax.set_xlabel("Amplitude (a.u.)")
    idx = Path(index).stem
    ax.set_title(
        f"{idx} — absolute amplitudes ({n} tidal factors)",
        fontsize=11,
    )
    ax.grid(True, axis="x", ls=":", lw=0.4, alpha=0.6)
    fig.tight_layout()
    plt.show()


def load_lte_results_time_forcing(target_dir: str) -> tuple[np.ndarray, np.ndarray]:
    """Load time (column 1) and forcing manifold (column 4) from lte_results.csv."""
    path = os.path.join(target_dir, "lte_results.csv")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"'lte_results.csv' not found in:\n  {target_dir}")

    try:
        data = np.loadtxt(path, delimiter=",", usecols=(0, 3), ndmin=2)
    except Exception as exc:
        raise ValueError(f"Failed to load lte_results.csv in:\n  {target_dir}\n\n{exc}") from exc

    if data.size == 0:
        raise ValueError(f"'lte_results.csv' is empty in:\n  {target_dir}")

    return data[:, 0], data[:, 1]


def plot_lpap_scatter_all(
    subdirs_info: list,
    start_year: int = 1950,
    n_years: int = 50,
) -> None:
    """
    Plot composed sinusoids for *all* datasets (one per subdir) as scatter
    dots (no connecting lines), so multiple overlapping series can be
    compared visually.

    Parameters
    ----------
    subdirs_info : list of (subdir_path: str, label: str)
    start_year   : first year of the time axis
    n_years      : length of the time axis in years
    """
    t_days, dates = _date_array(start_year, n_years, dt_days=1.0)

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.axhline(0, color="k", lw=0.4)

    n_total = max(len(subdirs_info), 1)
    colors = plt.cm.tab20(np.linspace(0, 1, n_total))

    loaded = 0
    for i, (subdir, label) in enumerate(subdirs_info):
        try:
            lpap = load_lpap(subdir)
            y = compose_sinusoids(lpap, t_days)
            ax.plot(
                dates, y,
                marker=".", markersize=1, linestyle="none",
                color=colors[i % n_total], alpha=0.6,
                label=label,
            )
            loaded += 1
        except Exception:
            pass

    ax.set_title(
        f"All datasets — composed tidal sinusoids "
        f"({loaded} of {len(subdirs_info)} loaded)\n"
        f"{start_year} – {start_year + n_years}",
        fontsize=11,
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Amplitude (a.u.)")
    ax.xaxis.set_major_locator(mdates.YearLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    plt.xticks(rotation=45, ha="right")
    ax.grid(True, ls=":", lw=0.4, alpha=0.6)
    fig.tight_layout()
    plt.show()


def plot_lte_results_forcing_all(subdirs_info: list) -> None:
    """
    Plot lte_results.csv column 1 (time) vs column 4 (forcing manifold)
    for all datasets, using very thin dash markers without connecting lines.
    """
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.axhline(0, color="k", lw=0.4)

    n_total = max(len(subdirs_info), 1)
    colors = plt.cm.tab20(np.linspace(0, 1, n_total))

    loaded = 0
    loaded_series: list[tuple[np.ndarray, np.ndarray]] = []
    for i, (subdir, _label) in enumerate(subdirs_info):
        try:
            times, forcing = load_lte_results_time_forcing(subdir)
            ax.plot(
                times, forcing,
                marker="_", markersize=2.0, markeredgewidth=0.25,
                linestyle="none",
                color=colors[i % n_total], alpha=0.7,
            )
            loaded_series.append((times, forcing))
            loaded += 1
        except Exception:
            pass

    if loaded == 0:
        raise ValueError("No valid lte_results.csv datasets found in any subdirectory.")

    avg_times = loaded_series[0][0]
    forcing_stack = np.full((loaded, len(avg_times)), np.nan)
    for i, (times, forcing) in enumerate(loaded_series):
        if len(times) == len(avg_times) and np.allclose(times, avg_times, rtol=0.0, atol=1e-9):
            forcing_stack[i] = forcing
        else:
            forcing_stack[i] = np.interp(avg_times, times, forcing, left=np.nan, right=np.nan)

    mean_forcing = np.nanmean(forcing_stack, axis=0)
    valid = ~np.isnan(mean_forcing)
    mean_path = Path(__file__).resolve().parent / "mean_forcing.dat"
    try:
        np.savetxt(mean_path, np.column_stack((avg_times[valid], mean_forcing[valid])), fmt="%.15e")
    except Exception as exc:
        raise RuntimeError(f"Failed to write mean forcing file:\n  {mean_path}\n\n{exc}") from exc

    ax.plot(avg_times[valid], mean_forcing[valid], "-", color="black", linewidth=0.6, zorder=10)

    ax.set_title(
        f"All datasets — forcing manifold from lte_results.csv "
        f"({loaded} of {len(subdirs_info)} loaded)",
        fontsize=11,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("Forcing manifold")
    ax.grid(True, ls=":", lw=0.4, alpha=0.6)
    fig.tight_layout()
    plt.show()


def plot_lpap_amplitudes_all(subdirs_info: list) -> None:
    """
    Plot tidal amplitude spectra for *all* datasets as outlined (no fill)
    horizontal bar charts overlaid on one axes, so amplitude scatter across
    datasets can be visualised.

    Parameters
    ----------
    subdirs_info : list of (subdir_path: str, label: str)
    """
    all_data: list = []
    for subdir, label in subdirs_info:
        try:
            lpap = load_lpap(subdir)
            periods    = [float(e[0]) for e in lpap]
            amplitudes = [abs(float(e[1])) for e in lpap]
            all_data.append((periods, amplitudes, label))
        except Exception:
            pass

    if not all_data:
        raise ValueError("No valid datasets found in any subdirectory.")

    # Build a union of all unique periods and map each to a y-position.
    all_periods_set: set = set()
    for periods, _, _ in all_data:
        all_periods_set.update(periods)
    all_periods = sorted(all_periods_set)

    n = len(all_periods)
    y_pos_map = {p: i for i, p in enumerate(all_periods)}
    y_pos     = np.arange(n)
    labels    = [f"{p:.5g} d" for p in all_periods]

    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * n)))

    n_total = max(len(all_data), 1)
    colors  = plt.cm.tab20(np.linspace(0, 1, n_total))

    for i, (periods, amplitudes, label) in enumerate(all_data):
        ys = [y_pos_map[p] for p in periods]
        ax.barh(
            ys, amplitudes,
            align="center",
            fill=False,
            edgecolor=colors[i % n_total],
            linewidth=0.8,
            label=label,
            alpha=0.8,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()

    def _sig5(x, _pos):
        if x == 0:
            return "0"
        from math import log10, floor
        digits = 5 - 1 - floor(log10(abs(x)))
        digits = max(0, digits)
        return f"{x:.{digits}f}"

    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(_sig5))

    ax.set_xlabel("Amplitude (a.u.)")
    ax.set_title(
        f"All datasets — tidal amplitude spectrum ({len(all_data)} datasets)",
        fontsize=11,
    )
    ax.grid(True, axis="x", ls=":", lw=0.4, alpha=0.6)
    fig.tight_layout()
    plt.show()


#######################################################################################

IS_WINDOWS = sys.platform == "win32"

# Optional PNG scaling
try:
    from PIL import Image, ImageTk  # type: ignore
    PIL_OK = True
except Exception:
    PIL_OK = False

# Optional YAML (preferred, since file is ID.yml)
try:
    import yaml  # type: ignore
    YAML_OK = True
except Exception:
    YAML_OK = False


# ---------------- CONFIG ----------------

PYTHON = sys.executable            # use the same interpreter running this script


# ---------------- HELPERS ----------------

def resolve_lt_cmd(run_dir: Path, use_json: bool = False) -> list[str]:
    """Return the lt command as an argument list."""
    json_flag = ["-j"] if use_json else []
    parent = run_dir.parent
    for name in ("lt.exe", "lt"):
        lt_path = parent / name
        if lt_path.exists():
            return [str(lt_path)] + json_flag
    raise FileNotFoundError(f"No lt or lt.exe found in {parent}")


def _find_terminal_linux() -> list[str] | None:
    """Return a terminal emulator prefix that can run a command, or None."""
    candidates = [
        (["xterm", "-e"],),
        (["gnome-terminal", "--"],),
        (["konsole", "-e"],),
        (["xfce4-terminal", "-e"],),
        (["lxterminal", "-e"],),
        (["mate-terminal", "-e"],),
    ]
    for (args,) in candidates:
        if shutil.which(args[0]):
            return args
    return None


def _open_terminal(cmd: list[str], cwd: Path, env: dict) -> None:
    """Open a new terminal window running *cmd* (list), keep it open after exit."""
    if IS_WINDOWS:
        subprocess.Popen(
            ["cmd.exe", "/k"] + cmd,
            cwd=str(cwd),
            env=env,
            creationflags=subprocess.CREATE_NEW_CONSOLE,
        )
        return

    title = str(cwd)
    set_title = f'echo -ne "\\033]0;{title}\\007"; '
    
    term = _find_terminal_linux()
    if term is None:
        # No GUI terminal found — run in background without a new window
        bare_cmd = " ".join(shlex.quote(a) for a in cmd)
        subprocess.Popen(
            ["bash", "-c", "ulimit -s unlimited && " + bare_cmd],
            cwd=str(cwd), env=env,
        )
        return

    term_name = term[0]
    # gnome-terminal uses -- to separate its args from the command
    shell_cmd = set_title + "ulimit -s unlimited && " + " ".join(shlex.quote(a) for a in cmd) + "; exec bash"
    subprocess.Popen(
        term + ["bash", "-c", shell_cmd],
        cwd=str(cwd), env=env,
    )


def load_id_map(id_file: Path) -> dict[str, dict]:
    """
    Load mapping of Index -> record from ID.yml (YAML) or JSON.
    Expected shapes supported:
      1) dict: { "14": {"Name": "...", "Country": "..."}, ... }
      2) list: [ {"Index":"14","Name":"...","Country":"..."}, ... ]
    """
    if not id_file.exists():
        return {}

    raw_text = id_file.read_text(encoding="utf-8", errors="replace")

    data = None

    # Prefer YAML if available; otherwise try JSON; finally attempt a small YAML subset.
    if YAML_OK:
        try:
            data = yaml.safe_load(raw_text)
        except Exception:
            data = None

    if data is None:
        try:
            data = json.loads(raw_text)
        except Exception:
            data = None

    # Tiny fallback for simple YAML maps if PyYAML isn't installed
    if data is None:
        # Very small subset: top-level keys with indented Name/Country:
        # 14:
        #   Name: Foo
        #   Country: Bar
        out: dict[str, dict] = {}
        cur_key: str | None = None
        for line in raw_text.splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if not line.startswith((" ", "\t")) and line.rstrip().endswith(":"):
                cur_key = line.rstrip()[:-1].strip().strip('"').strip("'")
                out[cur_key] = {}
                continue
            if cur_key and ":" in line:
                k, v = line.split(":", 1)
                k = k.strip()
                v = v.strip().strip('"').strip("'")
                out[cur_key][k] = v
        return out

    # Normalize to dict[str, dict]
    out: dict[str, dict] = {}

    if isinstance(data, dict):
        for k, v in data.items():
            out[str(k)] = v if isinstance(v, dict) else {"Value": v}
        return out

    if isinstance(data, list):
        for item in data:
            if isinstance(item, dict):
                key = item.get("Index") or item.get("ID") or item.get("index") or item.get("id")
                if key is not None:
                    out[str(key)] = item
        return out

    return {}


# ---------------- GUI ----------------

class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()

        self.title("LTE Runner (runs in selected index dir)")
        self.geometry("1100x780")

        self.root_dir: Path = Path.cwd().resolve()
        self.selected_index: str | None = None
        self.interval_begin_var = tk.StringVar(value="1940")
        self.interval_end_var = tk.StringVar(value="1970")
        self.metric_var = tk.StringVar(value="CC")  # internal values: "CC" or "DTW"

        # Load ID.yml from the same directory as this GUI script
        script_dir = Path(__file__).resolve().parent
        self.id_file = script_dir / "ID.yml"
        self.id_map: dict[str, dict] = load_id_map(self.id_file)

        self._image_ref = None
        self._pil_image_ref = None

        # README tooltip state
        self._tooltip_win: tk.Toplevel | None = None
        self._tooltip_after_id: str | None = None
        self._tooltip_item: str | None = None

        # JSON loading checkbox state
        self.use_json_var = tk.BooleanVar(value=True)  # Default to JSON mode

        self.exclude_var = tk.BooleanVar(value=True)  # Default to true
        self.filter_var = tk.BooleanVar(value=True)  # Default to true
        self.trend_var = tk.BooleanVar(value=True)  # Default to true
        self.test_only_var = tk.BooleanVar(value=False)  # Default to false
        self.sim_var = tk.BooleanVar(value=False)  # Use lte_results.csv col 2 as calibration target

        self._build_ui()
        self._set_root(self.root_dir)
        self._build_menu()
        
    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _build_menu(self):
        menubar = tk.Menu(self)

        # ── File menu ──────────────────────────────────────────────────
        file_menu = tk.Menu(menubar, tearoff=False)
        file_menu.add_command(
            label="Restore from github",
            # rm  lt.exe.p lt.exe.*.dat.p lt.exe.resp
            # git restore lt.exe.p lt.exe.*.dat.p lt.exe.resp
            command=self._cmd_restore_from_github,
        )
        file_menu.add_command(
            label="Copy from named index",
            # rm  lt.exe.p lt.exe.*.dat.p lt.exe.resp
            # cp ../INDEX_DIR/lt.exe.p ../INDEX_DIR/lt.exe.resp .
            command=self._cmd_copy_from_index,
        )
        file_menu.add_command(
            label="Edit RESP file",
            # gedit lt.exe.resp
            command=self._cmd_edit_resp,
        )
        file_menu.add_command(
            label="Edit JSON file",
            # gedit lt.exe.p
            command=self._cmd_edit_json,
        )
        file_menu.add_command(
            label="Remove JSON file",
            # gedit lt.exe.*.dat.p
            command=self._cmd_remove_json_dat,
        )
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        # ── Analyze menu ───────────────────────────────────────────────
        analyze_menu = tk.Menu(menubar, tearoff=False)
        analyze_menu.add_command(
            label="Plot Tidal Periodicities (lpap)…",
            command=self._cmd_plot_lpap,
        )
        analyze_menu.add_command(
            label="Plot Tidal Amplitude Spectrum (lpap)…",
            command=self._cmd_plot_lpap_amplitudes,
        )
        analyze_menu.add_separator()
        analyze_menu.add_command(
            label="Plot Tidal Periodicities — all subdirs (scatter)…",
            command=self._cmd_plot_lpap_all,
        )
        analyze_menu.add_command(
            label="Plot Tidal Amplitude Spectrum — all subdirs…",
            command=self._cmd_plot_lpap_amplitudes_all,
        )
        analyze_menu.add_command(
            label="Plot Forcing Manifold — all subdirs…",
            command=self._cmd_plot_lte_results_forcing_all,
        )
        menubar.add_cascade(label="Analyze", menu=analyze_menu)

        self.config(menu=menubar)


    # ------------------------------------------------------------------
    # Menu Callbacks
    # ------------------------------------------------------------------

    def _cmd_restore_from_github(self):
        """Menu command: File → Restore from github."""
        try:
            run_dir = self._run_dir()
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return
        if IS_WINDOWS:
            cmd = "del /F /Q lt.exe.p lt.exe.*.dat.p lt.exe.resp & git restore lt.exe.p lt.exe.*.dat.p lt.exe.resp"
            subprocess.Popen(["cmd.exe", "/c", cmd], cwd=str(run_dir))
        else:
            cmd = "rm -f lt.exe.p lt.exe.*.dat.p lt.exe.resp && git restore lt.exe.p lt.exe.*.dat.p lt.exe.resp"
            subprocess.Popen(["bash", "-c", cmd], cwd=str(run_dir))

    def _cmd_copy_from_index(self):
        """Menu command: File → Copy from named index."""
        try:
            run_dir = self._run_dir()
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return
        index_dir = simpledialog.askstring(
            "Copy from named index",
            "Enter source index directory name:",
            parent=self,
        )
        if not index_dir:
            return
        if IS_WINDOWS:
            src = index_dir.replace("/", "\\")
            cmd = (
                f"del /F /Q lt.exe.p lt.exe.*.dat.p lt.exe.resp"
                f" & copy ..\\{src}\\lt.exe.p . & copy ..\\{src}\\lt.exe.resp ."
            )
            subprocess.Popen(["cmd.exe", "/c", cmd], cwd=str(run_dir))
        else:
            src = shlex.quote(index_dir)
            cmd = f"rm -f lt.exe.p lt.exe.*.dat.p lt.exe.resp && cp ../{src}/lt.exe.p ../{src}/lt.exe.resp ."
            subprocess.Popen(["bash", "-c", cmd], cwd=str(run_dir))

    def _cmd_edit_resp(self):
        """Menu command: File → Edit RESP file."""
        try:
            run_dir = self._run_dir()
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return
        if IS_WINDOWS:
            subprocess.Popen(["notepad", "lt.exe.resp"], cwd=str(run_dir))
        else:
            subprocess.Popen(["gedit", "lt.exe.resp"], cwd=str(run_dir))

    def _cmd_edit_json(self):
        """Menu command: File → Edit JSON file."""
        try:
            run_dir = self._run_dir()
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return
        if IS_WINDOWS:
            subprocess.Popen(["notepad", "lt.exe.p"], cwd=str(run_dir))
        else:
            subprocess.Popen(["gedit", "lt.exe.p"], cwd=str(run_dir))

    def _cmd_remove_json_dat(self):
        """Menu command: File → Remove JSON file (edit lt.exe.*.dat.p)."""
        try:
            run_dir = self._run_dir()
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return
        if IS_WINDOWS:
            files = list(run_dir.glob("lt.exe.*.dat.p"))
            if not files:
                messagebox.showinfo("No files", "No lt.exe.*.dat.p files found.")
                return
            for f in files:
                subprocess.Popen(["notepad", str(f)])
        else:
            subprocess.Popen(["bash", "-c", "gedit lt.exe.*.dat.p"], cwd=str(run_dir))

    def _cmd_plot_lpap(self):
        """Menu command: Analyze → Plot Tidal Periodicities (lpap)."""
        try:
            target_dir = str(self._run_dir())
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return

        try:
            lpap = load_lpap(target_dir)
        except (FileNotFoundError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            messagebox.showerror("Error loading lt.exe.p", str(exc))
            return

        try:
            plot_lpap(lpap, target_dir, start_year=1950, n_years=50)
        except (ValueError, TypeError, RuntimeError, OverflowError) as exc:
            messagebox.showerror("Plot error", str(exc))
        

    def _cmd_plot_lpap_amplitudes(self):
        """Menu command: Analyze → Plot Tidal Amplitude Spectrum (lpap)."""
        try:
            target_dir = str(self._run_dir())
        except RuntimeError as exc:
            messagebox.showerror("No selection", str(exc))
            return

        try:
            lpap = load_lpap(target_dir)
        except (FileNotFoundError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            messagebox.showerror("Error loading lt.exe.p", str(exc))
            return

        try:
            plot_lpap_amplitudes(lpap, target_dir)
        except (ValueError, TypeError, RuntimeError, OverflowError) as exc:
            messagebox.showerror("Plot error", str(exc))

    # ------------------------------------------------------------------
    # All-subdirs helpers
    # ------------------------------------------------------------------

    def _get_all_subdirs_info(self, required_file: str = "lt.exe.p") -> list:
        """
        Return [(subdir_path, label), …] for every subdirectory under
        ``self.root_dir`` that contains *required_file*.

        Skips the standard non-data directories (locs, scripts, rlr_data).
        """
        SKIP = {"locs", "scripts", "rlr_data"}
        result = []
        try:
            for p in sorted(self.root_dir.iterdir()):
                if p.is_dir() and p.name not in SKIP and not p.name.startswith("qbo"):
                    if (p / required_file).exists():
                        result.append((str(p), p.name))
        except Exception as exc:
            raise RuntimeError(str(exc)) from exc
        return result

    def _cmd_plot_lpap_all(self):
        """Menu command: Analyze → Plot Tidal Periodicities — all subdirs (scatter)."""
        try:
            subdirs_info = self._get_all_subdirs_info()
        except RuntimeError as exc:
            messagebox.showerror("Error", str(exc))
            return

        if not subdirs_info:
            messagebox.showinfo(
                "No data",
                "No subdirectories with lt.exe.p were found under the current root."
            )
            return

        try:
            plot_lpap_scatter_all(subdirs_info, start_year=1950, n_years=50)
        except (ValueError, TypeError, RuntimeError, OverflowError) as exc:
            messagebox.showerror("Plot error", str(exc))

    def _cmd_plot_lpap_amplitudes_all(self):
        """Menu command: Analyze → Plot Tidal Amplitude Spectrum — all subdirs."""
        try:
            subdirs_info = self._get_all_subdirs_info()
        except RuntimeError as exc:
            messagebox.showerror("Error", str(exc))
            return

        if not subdirs_info:
            messagebox.showinfo(
                "No data",
                "No subdirectories with lt.exe.p were found under the current root."
            )
            return

        try:
            plot_lpap_amplitudes_all(subdirs_info)
        except (ValueError, TypeError, RuntimeError, OverflowError) as exc:
            messagebox.showerror("Plot error", str(exc))

    def _cmd_plot_lte_results_forcing_all(self):
        """Menu command: Analyze → Plot Forcing Manifold — all subdirs."""
        try:
            subdirs_info = self._get_all_subdirs_info(required_file="lte_results.csv")
        except RuntimeError as exc:
            messagebox.showerror("Error", str(exc))
            return

        if not subdirs_info:
            messagebox.showinfo(
                "No data",
                "No subdirectories with lte_results.csv were found under the current root."
            )
            return

        try:
            plot_lte_results_forcing_all(subdirs_info)
        except (ValueError, TypeError, RuntimeError, OverflowError) as exc:
            messagebox.showerror("Plot error", str(exc))

    # ---------- UI ----------

    def _build_ui(self) -> None:
        top = ttk.Frame(self, padding=8)
        top.pack(fill="x")

        ttk.Button(top, text="Choose Root…", command=self.choose_root).pack(side="left")
        ttk.Button(top, text="Refresh", command=self.refresh_list).pack(side="left", padx=8)

        ttk.Label(top, text="Root:").pack(side="left", padx=8)
        self.root_var = tk.StringVar()
        ttk.Entry(top, textvariable=self.root_var, width=90).pack(side="left", fill="x", expand=True)

        main = ttk.PanedWindow(self, orient="horizontal")
        main.pack(fill="both", expand=True, padx=8, pady=8)

        left = ttk.Frame(main, padding=6)
        right = ttk.Frame(main, padding=6)
        main.add(left, weight=1)
        main.add(right, weight=10)

        ttk.Label(left, text="Index directories (select one):").pack(anchor="w")
        self.dir_list = tk.Listbox(left, height=28)
        self.dir_list.pack(fill="both", expand=True)
        self.dir_list.bind("<<ListboxSelect>>", lambda e: self._on_select())
        self.dir_list.bind("<Motion>", self._on_list_motion)
        self.dir_list.bind("<Leave>", self._on_list_leave)

        sel_frame = ttk.LabelFrame(right, text="Selected directory + ID.yml info", padding=8)
        sel_frame.pack(fill="x")

        self.sel_var = tk.StringVar(value="(none)")
        ttk.Entry(sel_frame, textvariable=self.sel_var, state="readonly").pack(fill="x")

        info_grid = ttk.Frame(sel_frame)
        info_grid.pack(fill="x", pady=(8, 0))

        ttk.Label(info_grid, text="Name:").grid(row=0, column=0, sticky="w")
        self.name_var = tk.StringVar(value="")
        ttk.Entry(info_grid, textvariable=self.name_var, state="readonly").grid(row=0, column=1, sticky="we", padx=(8, 0))

        ttk.Label(info_grid, text="Country:").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.country_var = tk.StringVar(value="")
        ttk.Entry(info_grid, textvariable=self.country_var, state="readonly").grid(row=1, column=1, sticky="we", padx=(8, 0), pady=(6, 0))

        info_grid.columnconfigure(1, weight=1)

        run_frame = ttk.LabelFrame(right, text="Run (interactive console)", padding=8)
        run_frame.pack(fill="x", pady=10)

        ttk.Label(run_frame, text="TIMEOUT (seconds):").grid(row=0, column=0, sticky="w")
        self.time_var = tk.DoubleVar(value=1000.0)
        ttk.Scale(run_frame, from_=0.0, to=36000.0, variable=self.time_var).grid(
            row=0, column=1, sticky="we", padx=8
        )
        run_frame.columnconfigure(1, weight=1)

        self.time_lbl = ttk.Label(run_frame, text=f"{self.time_var.get():,.3f}")
        self.time_lbl.grid(row=0, column=2, sticky="e")
        self.time_var.trace_add("write", lambda *_: self.time_lbl.config(text=f"{self.time_var.get():,.3f}"))

        btns = ttk.Frame(run_frame)
        btns.grid(row=1, column=0, columnspan=3, pady=10, sticky="w")

        # ttk.Button(btns, text="Run lt", command=self.run_lt).pack(side="left")
        # ttk.Button(btns, text="Run plot", command=self.run_plot).pack(side="left", padx=8)
        # ttk.Button(btns, text="Refresh PNG", command=self.show_png).pack(side="left", padx=8)

        btns = ttk.Frame(run_frame)
        btns.grid(row=1, column=0, columnspan=3, pady=10, sticky="we")
        btns.columnconfigure(0, weight=1)

        # Left side: buttons
        left_btns = ttk.Frame(btns)
        left_btns.grid(row=0, column=0, sticky="w")

        ttk.Button(left_btns, text="Run lt", command=self.run_lt).pack(side="left")
        ttk.Checkbutton(left_btns, text="JSON", variable=self.use_json_var).pack(side="left", padx=(8, 0))
        ttk.Button(left_btns, text="Run plot", command=self.run_plot).pack(side="left", padx=8)
        ttk.Button(left_btns, text="Refresh PNG", command=self.show_png).pack(side="left", padx=8)

        # Right side: interval begin/end editors (LAST on the row)
        # right_fields = ttk.Frame(btns)
        # right_fields.grid(row=0, column=1, sticky="e")

        # ttk.Label(right_fields, text="Interval:").pack(side="left", padx=(8, 6))
        # ttk.Entry(right_fields, textvariable=self.interval_begin_var, width=6).pack(side="left")
        # ttk.Label(right_fields, text="to").pack(side="left", padx=6)
        # ttk.Entry(right_fields, textvariable=self.interval_end_var, width=6).pack(side="left")

        right_fields = ttk.Frame(btns)
        right_fields.grid(row=0, column=1, sticky="e")

        # Metric radios (right before Interval boxes)
        ttk.Label(right_fields, text="Metric:").pack(side="left", padx=(8, 6))

        ttk.Radiobutton(
            right_fields, text="CC", variable=self.metric_var, value="CC"
        ).pack(side="left")

        ttk.Radiobutton(
            right_fields, text="DTW", variable=self.metric_var, value="DTW"
        ).pack(side="left", padx=(6, 0))

        ttk.Radiobutton(
            right_fields, text="H", variable=self.metric_var, value="HOYER"
        ).pack(side="left", padx=(6, 0))


        # Interval editors (LAST on the row)
#        ttk.Label(right_fields, text="  Test Interval:").pack(side="left", padx=(12, 6))
        ttk.Label(right_fields, text="  CV:").pack(side="left", padx=(12, 6))
        ttk.Entry(right_fields, textvariable=self.interval_begin_var, width=6).pack(side="left")
        ttk.Label(right_fields, text="to").pack(side="left", padx=6)
        ttk.Entry(right_fields, textvariable=self.interval_end_var, width=6).pack(side="left")
        ttk.Checkbutton(right_fields, text="exclude", variable=self.exclude_var).pack(side="left")
        ttk.Checkbutton(right_fields, text="filter", variable=self.filter_var).pack(side="left")
        ttk.Checkbutton(right_fields, text="trend", variable=self.trend_var).pack(side="left")
        ttk.Checkbutton(right_fields, text="test", variable=self.test_only_var).pack(side="left")
        ttk.Checkbutton(right_fields, text="sim", variable=self.sim_var).pack(side="left")


        img_frame = ttk.LabelFrame(right, text="PNG preview (from selected dir)", padding=8)
        img_frame.pack(fill="both", expand=True)

        self.canvas = tk.Canvas(img_frame, background="black")
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", lambda e: self._redraw_image())

    # ---------- Root / selection ----------

    def choose_root(self) -> None:
        path = filedialog.askdirectory(title="Choose root folder (contains index dirs)")
        if path:
            self._set_root(Path(path).resolve())

    def _set_root(self, root: Path) -> None:
        self.root_dir = root
        self.root_var.set(str(root))
        self.selected_index = None
        self.sel_var.set("(none)")
        self.name_var.set("")
        self.country_var.set("")
        self._clear_image()
        self.refresh_list()

    def _on_list_motion(self, event: tk.Event) -> None:
        idx = self.dir_list.nearest(event.y)
        if idx < 0 or idx >= self.dir_list.size():
            self._hide_readme_tooltip()
            return
        item = self.dir_list.get(idx)
        if item == self._tooltip_item:
            return
        self._hide_readme_tooltip()
        self._tooltip_item = item
        readme_path = self.root_dir / item / "README.md"
        if not readme_path.is_file():
            return
        try:
            text = readme_path.read_text(encoding="utf-8").strip()
        except Exception:
            return
        if not text:
            return
        # Schedule the tooltip to appear after a short delay
        self._tooltip_after_id = self.after(400, lambda: self._show_readme_tooltip(text, event))

    def _on_list_leave(self, event: tk.Event) -> None:
        self._hide_readme_tooltip()

    def _show_readme_tooltip(self, text: str, event: tk.Event) -> None:
        self._hide_readme_tooltip(cancel_only=True)
        win = tk.Toplevel(self)
        win.wm_overrideredirect(True)
        win.wm_attributes("-topmost", True)
        win.wm_geometry(f"+{event.x_root + 20}+{event.y_root + 10}")
        lbl = tk.Label(
            win,
            text=text,
            justify="left",
            background="#ffffe0",
            relief="solid",
            borderwidth=1,
            wraplength=420,
            padx=6,
            pady=4,
        )
        lbl.pack()
        self._tooltip_win = win

    def _hide_readme_tooltip(self, cancel_only: bool = False) -> None:
        if self._tooltip_after_id is not None:
            self.after_cancel(self._tooltip_after_id)
            self._tooltip_after_id = None
        if not cancel_only and self._tooltip_win is not None:
            self._tooltip_win.destroy()
            self._tooltip_win = None
            self._tooltip_item = None

    def refresh_list(self) -> None:
        self.dir_list.delete(0, tk.END)
        try:
            for p in sorted(self.root_dir.iterdir()):
                # if p.is_dir():
                if p.is_dir() and p.name not in {"locs", "scripts", "rlr_data"}:
                    self.dir_list.insert(tk.END, p.name)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _on_select(self) -> None:
        sel = self.dir_list.curselection()
        if not sel:
            self.selected_index = None
            self.sel_var.set("(none)")
            self.name_var.set("")
            self.country_var.set("")
            return

        self.selected_index = self.dir_list.get(sel[0])
        self.sel_var.set(str((self.root_dir / self.selected_index).resolve()))
        self._clear_image()
        self._update_id_fields()
        self.show_loc_png_for_selected()

    def _run_dir(self) -> Path:
        if not self.selected_index:
            raise RuntimeError("Select an index directory first.")
        return (self.root_dir / self.selected_index).resolve()

    def _update_id_fields(self) -> None:
        idx = (self.selected_index or "").strip().strip('"').strip("'")
        rec = self.id_map.get(idx)

        name = ""
        country = ""

        if isinstance(rec, dict):
            name = str(rec.get("Name", "") or rec.get("name", "") or "")
            country = str(rec.get("Country", "") or rec.get("country", "") or "")

        self.name_var.set(name)
        self.country_var.set(country)

    # ---------- Actions ----------

    def run_lt(self) -> None:
        try:
            run_dir = self._run_dir()
        except Exception as e:
            messagebox.showerror("No selection", str(e))
            return

        index = run_dir.name.strip().strip('"').strip("'")
        timeout = float(self.time_var.get())

        par = run_dir / f"lt.exe.{index}.dat.par"
        if par.exists():
            try:
                par.unlink()
            except Exception:
                pass

        try:
            use_json = self.use_json_var.get()
            lt_cmd = resolve_lt_cmd(run_dir, use_json)
        except Exception as e:
            messagebox.showerror("Missing lt", str(e))
            return
        b = self.interval_begin_var.get().strip()
        e = self.interval_end_var.get().strip()
        if not b or not e:
            messagebox.showerror("Bad interval", "Interval begin/end must be non-empty.")
            return

        env = os.environ.copy()
        env["METRIC"] = self.metric_var.get().strip().upper()
        env["TIMEOUT"] = f"{timeout:.6f}"
        env["TRAIN_START"] = b
        env["TRAIN_END"] = e
        env["CLIMATE_INDEX"] = f"{index}.dat"
        if self.sim_var.get():
            csv_path = run_dir / "lte_results.csv"
            if not csv_path.exists():
                messagebox.showerror(
                    "sim mode",
                    f"lte_results.csv not found in:\n  {run_dir}\n\n"
                    "Run the model once first to generate it."
                )
                return
            sim_dat = run_dir / "lte_results_model.dat"
            try:
                with csv_path.open() as fh, sim_dat.open("w") as out:
                    for line in fh:
                        line = line.strip()
                        if not line:
                            continue
                        parts = [p.strip() for p in line.split(",")]
                        if len(parts) >= 2:
                            out.write(f"{parts[0]}  {parts[1]}\n")
            except Exception as exc:
                messagebox.showerror("sim mode", f"Failed to create sim data file:\n{exc}")
                return
            env["CLIMATE_INDEX"] = sim_dat.name
        env["IDATE"] = "1920.9"
        env["EXCLUDE"] = "true" if self.exclude_var.get() else "false"
        env["TREND"] = "true" if self.trend_var.get() else "false"
        env["F9"] = "1" if self.filter_var.get() else "0"
        env["TEST_ONLY"] = "true" if self.test_only_var.get() else "false"

        _open_terminal(lt_cmd, run_dir, env)

    def run_plot(self) -> None:
        try:
            run_dir = self._run_dir()
        except Exception as e:
            messagebox.showerror("No selection", str(e))
            return

        index = run_dir.name.strip().strip('"').strip("'")

        plot_path = (run_dir.parent / "plot.py").resolve()
        if not plot_path.exists():
            messagebox.showerror(
                "Missing plot.py",
                f"Could not find plot.py at:\n{plot_path}"
            )
            return

        begin = self.interval_begin_var.get().strip()
        end   = self.interval_end_var.get().strip()

        if IS_WINDOWS:
            cmd = [PYTHON, str(plot_path), index, "Feb2026", begin, end, "0"]
            _open_terminal(cmd, run_dir, os.environ.copy())
            return

        # On Linux: execute plot.py in a subprocess so that matplotlib starts
        # fresh with the Agg backend (the GUI process already loaded TkAgg on
        # the main thread, so running plot.py in-process via runpy would trigger
        # "Starting a Matplotlib GUI outside of the main thread" warnings).
        # The subprocess is launched from a daemon thread so the GUI stays
        # responsive; the PNG preview is auto-refreshed when the script finishes.
        def _run() -> None:
            env = os.environ.copy()
            env["MPLBACKEND"] = "Agg"  # non-interactive; prevents any window
            cmd = [PYTHON, str(plot_path), index, "Feb2026", begin, end, "0"]
            try:
                result = subprocess.run(
                    cmd,
                    cwd=str(run_dir),
                    env=env,
                    capture_output=True,
                    text=True,
                )
                if result.returncode != 0:
                    err = (result.stderr or result.stdout).strip()
                    self.after(0, lambda err=err: messagebox.showerror("Plot error", err))
            except Exception as exc:
                self.after(0, lambda exc=exc: messagebox.showerror("Plot error", str(exc)))

            self.after(0, self.show_png)

        threading.Thread(target=_run, daemon=True).start()

    # ---------- PNG ----------

    def show_png(self) -> None:
        try:
            run_dir = self._run_dir()
        except Exception as e:
            messagebox.showerror("No selection", str(e))
            return

        pngs = list(run_dir.glob("*.png"))
        if not pngs:
            messagebox.showinfo("No PNG", "No PNG found in selected directory.")
            return

        pngs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        self._load_image(pngs[0])

    def _clear_image(self) -> None:
        self.canvas.delete("all")
        self._image_ref = None
        self._pil_image_ref = None

    def _load_image(self, path: Path) -> None:
        self._clear_image()
        try:
            if PIL_OK:
                img = Image.open(path)
                self._pil_image_ref = img
                self._redraw_image()
            else:
                self._image_ref = tk.PhotoImage(file=str(path))
                self.canvas.create_image(0, 0, image=self._image_ref, anchor="nw")
        except Exception as e:
            messagebox.showerror("Image error", str(e))

    def _redraw_image(self) -> None:
        if not PIL_OK or self._pil_image_ref is None:
            return

        self.canvas.delete("all")
        cw = max(1, self.canvas.winfo_width())
        ch = max(1, self.canvas.winfo_height())

        img = self._pil_image_ref.copy()
        iw, ih = img.size
        scale = min(cw / iw, ch / ih)
        img = img.resize((int(iw * scale), int(ih * scale)))

        self._image_ref = ImageTk.PhotoImage(img)
        self.canvas.create_image(cw // 2, ch // 2, image=self._image_ref, anchor="center")

    def show_loc_png_for_selected(self) -> None:
        """
        Auto-preview: show locs/<index>_loc.png for the selected index.
        Expected location relative to selected index dir:
          ../locs/<index>_loc.png
        """
        try:
            run_dir = self._run_dir()
        except Exception:
            return

        index = run_dir.name.strip().strip('"').strip("'")

        # Primary intended location: parent/locs/<index>_loc.png
        candidates = [
            run_dir.parent / "locs" / f"{index}_loc.png",
            run_dir / "locs" / f"{index}_loc.png",            # fallback
            self.root_dir / "locs" / f"{index}_loc.png",      # fallback
        ]

        for p in candidates:
            if p.exists():
                self._load_image(p)
                return

        # If not found, just clear the preview (don’t popup)
        self._clear_image()

def main() -> None:
    App().mainloop()


if __name__ == "__main__":
    main()
