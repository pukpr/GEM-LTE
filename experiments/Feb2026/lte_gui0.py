#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
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


# ---------------- CONFIG ----------------

PYTHON = "python3"                 # Windows launcher
PLOT_REL = r"..\plot.py"      # relative to INDEX directory


# ---------------- HELPERS ----------------

def resolve_lt_cmd(run_dir: Path) -> str:
    if (run_dir / "lt.exe").exists():
        return r".\lt.exe"
    if (run_dir / "lt").exists():
        return r".\lt"
    raise FileNotFoundError(f"No lt or lt.exe in {run_dir}")


# ---------------- GUI ----------------

class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()

        self.title("LTE Runner (runs in selected index dir)")
        self.geometry("1100x740")

        self.root_dir: Path = Path.cwd().resolve()
        self.selected_index: str | None = None

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
        main.add(right, weight=2)

        ttk.Label(left, text="Index directories (select one):").pack(anchor="w")
        self.dir_list = tk.Listbox(left, height=28)
        self.dir_list.pack(fill="both", expand=True)
        self.dir_list.bind("<<ListboxSelect>>", lambda e: self._on_select())

        sel_frame = ttk.LabelFrame(right, text="Selected directory", padding=8)
        sel_frame.pack(fill="x")

        self.sel_var = tk.StringVar(value="(none)")
        ttk.Entry(sel_frame, textvariable=self.sel_var, state="readonly").pack(fill="x")

        run_frame = ttk.LabelFrame(right, text="Run (interactive console)", padding=8)
        run_frame.pack(fill="x", pady=10)

        ttk.Label(run_frame, text="TIMEOUT (seconds):").grid(row=0, column=0, sticky="w")
        self.time_var = tk.DoubleVar(value=1000.0)
        ttk.Scale(run_frame, from_=0.0, to=10000.0, variable=self.time_var).grid(
            row=0, column=1, sticky="we", padx=8
        )
        run_frame.columnconfigure(1, weight=1)

        self.time_lbl = ttk.Label(run_frame, text=f"{self.time_var.get():,.3f}")
        self.time_lbl.grid(row=0, column=2, sticky="e")
        self.time_var.trace_add("write", lambda *_: self.time_lbl.config(text=f"{self.time_var.get():,.3f}"))

        btns = ttk.Frame(run_frame)
        btns.grid(row=1, column=0, columnspan=3, pady=10, sticky="w")

        ttk.Button(btns, text="Run lt", command=self.run_lt).pack(side="left")
        ttk.Button(btns, text="Run plot", command=self.run_plot).pack(side="left", padx=8)
        ttk.Button(btns, text="Refresh PNG", command=self.show_png).pack(side="left", padx=8)

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
        self._clear_image()
        self.refresh_list()

    def refresh_list(self) -> None:
        self.dir_list.delete(0, tk.END)
        try:
            for p in sorted(self.root_dir.iterdir()):
                if p.is_dir():
                    self.dir_list.insert(tk.END, p.name)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _on_select(self) -> None:
        sel = self.dir_list.curselection()
        if not sel:
            self.selected_index = None
            self.sel_var.set("(none)")
            return
        self.selected_index = self.dir_list.get(sel[0])
        self.sel_var.set(str((self.root_dir / self.selected_index).resolve()))
        self._clear_image()

    def _run_dir(self) -> Path:
        if not self.selected_index:
            raise RuntimeError("Select an index directory first.")
        return (self.root_dir / self.selected_index).resolve()

    # ---------- Actions ----------

    def run_lt(self) -> None:
        try:
            run_dir = self._run_dir()
        except Exception as e:
            messagebox.showerror("No selection", str(e))
            return

        index = run_dir.name
        timeout = float(self.time_var.get())

        # clean param file
        par = run_dir / f"lt.exe.{index}.dat.par"
        if par.exists():
            try:
                par.unlink()
            except Exception:
                pass

        try:
            lt_cmd = resolve_lt_cmd(run_dir)
        except Exception as e:
            messagebox.showerror("Missing lt", str(e))
            return

        env = os.environ.copy()
        env["METRIC"] = "cc"
        env["TIMEOUT"] = f"{timeout:.6f}"
        env["TRAIN_START"] = "1940"
        env["TRAIN_STOP"] = "1970"
        env["CLIMATE_INDEX"] = f"{index}.dat"
        env["IDATE"] = "1920.9"

        subprocess.Popen(
            ["cmd.exe", "/k", lt_cmd],
            cwd=str(run_dir),              # <<< THIS IS THE KEY
            env=env,
            creationflags=subprocess.CREATE_NEW_CONSOLE,
        )


    def run_plot(self) -> None:
        try:
            run_dir = self._run_dir()
        except Exception as e:
            messagebox.showerror("No selection", str(e))
            return

        index = run_dir.name.strip().strip('"').strip("'")  # sanitize

        # Resolve plot.py relative to run_dir
        plot_path = (run_dir / PLOT_REL).resolve()
        if not plot_path.exists():
            messagebox.showerror("Missing plot.py", f"Could not find plot.py at:\n{plot_path}")
            return

        args = [
            PYTHON, str(plot_path),
            index, "Feb2026", "1940", "1970", "0"
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


def main() -> None:
    App().mainloop()


if __name__ == "__main__":
    main()

