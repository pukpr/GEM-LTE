#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import json
import subprocess
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

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

PYTHON = "python3"                 # Windows launcher
PLOT_REL = r"..\plot.py"      # relative to INDEX directory


# ---------------- HELPERS ----------------

def resolve_lt_cmd(run_dir: Path) -> str:
    if (run_dir / "..\\lt.exe").exists():
        return r"..\\lt.exe"
    if (run_dir / "..\\lt").exists():
        return r"..\\lt"
    raise FileNotFoundError(f"No lt or lt.exe in {run_dir}")


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

        self._build_ui()
        self._set_root(self.root_dir)

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

        # Interval editors (LAST on the row)
        ttk.Label(right_fields, text="  Test Interval:").pack(side="left", padx=(12, 6))
        ttk.Entry(right_fields, textvariable=self.interval_begin_var, width=6).pack(side="left")
        ttk.Label(right_fields, text="to").pack(side="left", padx=6)
        ttk.Entry(right_fields, textvariable=self.interval_end_var, width=6).pack(side="left")


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
            lt_cmd = resolve_lt_cmd(run_dir)
            # lt_cmd = resolve_lt_cmd(run_dir.parent)
        except Exception as e:
            messagebox.showerror("Missing lt", str(e))
            return
        b = self.interval_begin_var.get().strip()
        e = self.interval_end_var.get().strip()
        if not b or not e:
            messagebox.showerror("Bad interval", "Interval begin/end must be non-empty.")
            return

        env = os.environ.copy()
        # env["METRIC"] = "cc"
        env["METRIC"] = self.metric_var.get().strip().upper()
        env["TIMEOUT"] = f"{timeout:.6f}"
        # env["TRAIN_START"] = "1940"
        # env["TRAIN_STOP"] = "1970"
        env["TRAIN_START"] = b # self.interval_begin_var.get().strip()
        env["TRAIN_STOP"] = e # self.interval_end_var.get().strip()
        env["CLIMATE_INDEX"] = f"{index}.dat"
        env["IDATE"] = "1920.9"

        subprocess.Popen(
            ["cmd.exe", "/k", lt_cmd],
            cwd=str(run_dir),
            env=env,
            creationflags=subprocess.CREATE_NEW_CONSOLE,
        )

    def run_plot(self) -> None:
        try:
            run_dir = self._run_dir()
        except Exception as e:
            messagebox.showerror("No selection", str(e))
            return

        index = run_dir.name.strip().strip('"').strip("'")

        plot_path = (run_dir / PLOT_REL).resolve()
        if not plot_path.exists():
            messagebox.showerror(
                "Missing plot.py",
                f"Could not find plot.py at:\n{plot_path}\n\nAdjust PLOT_REL at top of script."
            )
            return

        args = [
            PYTHON,
            str(plot_path),
            index,
            "Feb2026",
            # "1940",
            # "1970",
            self.interval_begin_var.get().strip(),
            self.interval_end_var.get().strip(),
            "0",
        ]

        subprocess.Popen(
            args,
            cwd=str(run_dir),
            creationflags=subprocess.CREATE_NEW_CONSOLE,
        )

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
