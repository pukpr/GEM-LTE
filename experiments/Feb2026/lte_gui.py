#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import sys
import json
import runpy
import shlex
import shutil
import subprocess
import threading
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import datetime


# ---------------------------------------------------------------------------
# Core computation
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


def plot_lpap(lpap: list, start_year: int = 1950, n_years: int = 50):
    """
    Plot the composed sum of sinusoids defined by *lpap* over *n_years*
    starting at *start_year*.
    """
    t_days, dates = _date_array(start_year, n_years, dt_days=1.0)
    y = compose_sinusoids(lpap, t_days)

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(dates, y, lw=0.8, color="steelblue")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_title(
        f"Tidal Periodicities (lpap) — composed sum of {len(lpap)} sinusoid(s)\n"
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
        
        # JSON loading checkbox state
        self.use_json_var = tk.BooleanVar(value=True)  # Default to JSON mode

        self.exclude_var = tk.BooleanVar(value=True)  # Default to true
        self.filter_var = tk.BooleanVar(value=True)  # Default to true
        self.trend_var = tk.BooleanVar(value=True)  # Default to true
        self.test_only_var = tk.BooleanVar(value=False)  # Default to false

        self._build_ui()
        self._set_root(self.root_dir)
        self._build_menu()
        # self._build_body()
        
    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _build_menu(self):
        menubar = tk.Menu(self)

        # ── File menu ──────────────────────────────────────────────────
        file_menu = tk.Menu(menubar, tearoff=False)
        file_menu.add_command(
            label="Set target directory…",
            command=self._choose_directory,
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
        menubar.add_cascade(label="Analyze", menu=analyze_menu)

        self.config(menu=menubar)

    def _build_body(self):
        pad = {"padx": 10, "pady": 6}
        tk.Label(self, text="Target directory:").grid(
            row=0, column=0, sticky="w", **pad
        )
        tk.Label(
            self,
            textvariable=self._target_dir,
            relief="sunken",
            width=42,
            anchor="w",
        ).grid(row=0, column=1, sticky="ew", **pad)
        tk.Button(
            self, text="Browse…", command=self._choose_directory
        ).grid(row=1, column=1, sticky="e", **pad)
        self.columnconfigure(1, weight=1)

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _choose_directory(self):
        directory = filedialog.askdirectory(title="Select target directory")
        if directory:
            self._target_dir.set(directory)

    def _cmd_plot_lpap(self):
        """Menu command: Analyze → Plot Tidal Periodicities (lpap)."""
        target_dir = self._target_dir.get()
        if target_dir == "(no directory selected)":
            # Prompt the user to pick a directory first
            directory = filedialog.askdirectory(
                title="Select directory containing lt.exe.p"
            )
            if not directory:
                return
            self._target_dir.set(directory)
            target_dir = directory

        try:
            lpap = load_lpap(target_dir)
        except (FileNotFoundError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            messagebox.showerror("Error loading lt.exe.p", str(exc))
            return

        try:
            plot_lpap(lpap, start_year=1950, n_years=50)
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

        # On Linux: execute plot.py in-process (no terminal window).
        # Runs in a background thread so the GUI stays responsive; auto-refreshes
        # the PNG preview when the script finishes.
        def _run() -> None:
            saved_argv = sys.argv[:]
            saved_cwd  = os.getcwd()
            saved_mpl  = os.environ.get("MPLBACKEND")

            sys.argv = [str(plot_path), index, "Feb2026", begin, end, "0"]
            os.chdir(str(run_dir))
            os.environ["MPLBACKEND"] = "Agg"  # non-interactive; prevents any window

            try:
                runpy.run_path(str(plot_path), run_name="__main__")
            except SystemExit:
                pass
            except Exception as exc:
                self.after(0, lambda exc=exc: messagebox.showerror("Plot error", str(exc)))
            finally:
                sys.argv = saved_argv
                os.chdir(saved_cwd)
                if saved_mpl is None:
                    os.environ.pop("MPLBACKEND", None)
                else:
                    os.environ["MPLBACKEND"] = saved_mpl
                # Release any matplotlib figures created by plot.py
                try:
                    import matplotlib.pyplot as _plt
                    _plt.close("all")
                except Exception:
                    pass

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
