#!/usr/bin/env python3
import os
import json
import math
from collections import defaultdict
from typing import Any, Dict, List, Tuple, Iterable, Optional

import matplotlib.pyplot as plt

# -----------------------
# CONFIG YOU MAY TWEAK
# -----------------------
TARGET_FILENAME = "ts.dat.p"     # expected JSON file name in each subdirectory
OUTPUT_PREFIX = "hist_page"      # output PNGs: hist_page_001.png, ...

DPI = 300                        # high DPI for zooming
MAX_COLS = 6
MAX_ROWS = 6
BINS = 30

# Keys used to detect paired amplitude/phase/frequency lists inside the same dict.
# Add variants you use in your JSON.
AMP_KEYS   = {"amp", "amps", "amplitude", "amplitudes"}
PHASE_KEYS = {"phase", "phases", "phi", "phasor_phase"}
FREQ_KEYS  = {"freq", "freqs", "frequency", "frequencies", "stable_freq", "stable_freqs"}

TWO_PI = 2.0 * math.pi
PI = math.pi


def is_float(x: Any) -> bool:
    return isinstance(x, float) and not isinstance(x, bool)


def is_float_list(x: Any) -> bool:
    return isinstance(x, list) and all(is_float(v) for v in x)


def wrap_0_2pi(x: float) -> float:
    # Wrap to [0, 2π)
    y = x % TWO_PI
    # % already yields [0, TWO_PI) for positive TWO_PI, but keep it explicit
    if y < 0:
        y += TWO_PI
    return y


def find_first_key(d: Dict[str, Any], keyset: set) -> Optional[str]:
    """Return the first matching key (in insertion order) from keyset, or None."""
    for k in d.keys():
        if k in keyset:
            return k
    return None


def correct_amp_phase_inplace(obj: Any) -> None:
    """
    Walk the object tree. Wherever we find a dict containing amplitude list + phase list
    (and optionally frequency list) of same length, apply the sign/phase correction:
      if amp[i] < 0: amp[i] = -amp[i], phase[i] = (phase[i] + π) mod 2π
      phase[i] always wrapped to [0, 2π)
    """
    if isinstance(obj, dict):
        # Recurse first (so nested dicts corrected too)
        for v in obj.values():
            correct_amp_phase_inplace(v)

        amp_k = find_first_key(obj, AMP_KEYS)
        ph_k = find_first_key(obj, PHASE_KEYS)
        if amp_k is None or ph_k is None:
            return

        amp = obj.get(amp_k)
        ph = obj.get(ph_k)
        if not (is_float_list(amp) and is_float_list(ph)):
            return
        if len(amp) != len(ph):
            return

        # Apply correction
        for i in range(len(amp)):
            if amp[i] < 0.0:
                amp[i] = -amp[i]
                ph[i] = ph[i] + PI
            ph[i] = wrap_0_2pi(ph[i])

        # write back (lists are mutated anyway, but keep it explicit)
        obj[amp_k] = amp
        obj[ph_k] = ph

        # No need to touch frequency list; it's used only for labeling later.

    elif isinstance(obj, list):
        for v in obj:
            correct_amp_phase_inplace(v)


def build_freq_label_map(obj: Any) -> Dict[Tuple[str, int], float]:
    """
    Build a map from (prefix_path + '.<amp_or_phase_key>', index) -> frequency_value
    whenever we find sibling stable frequency list in the same dict.

    Example: if at prefix 'model.modes' we have:
      {"amplitude":[...], "phase":[...], "stable_freqs":[...]}
    then we create entries mapping:
      ('model.modes.amplitude', i) -> stable_freqs[i]
      ('model.modes.phase', i)     -> stable_freqs[i]
    """
    freq_map: Dict[Tuple[str, int], float] = {}

    def walk(x: Any, prefix: str = ""):
        if isinstance(x, dict):
            # Identify siblings if present
            amp_k = find_first_key(x, AMP_KEYS)
            ph_k = find_first_key(x, PHASE_KEYS)
            fq_k = find_first_key(x, FREQ_KEYS)

            amp = x.get(amp_k) if amp_k else None
            ph = x.get(ph_k) if ph_k else None
            fq = x.get(fq_k) if fq_k else None

            if fq_k and is_float_list(fq):
                # If amplitudes present and same length, map for labeling
                if amp_k and is_float_list(amp) and len(amp) == len(fq):
                    base = f"{prefix}.{amp_k}" if prefix else amp_k
                    for i, fval in enumerate(fq):
                        freq_map[(base, i)] = fval
                # Same for phases
                if ph_k and is_float_list(ph) and len(ph) == len(fq):
                    base = f"{prefix}.{ph_k}" if prefix else ph_k
                    for i, fval in enumerate(fq):
                        freq_map[(base, i)] = fval

            # Recurse
            for k, v in x.items():
                k_str = str(k)
                new_prefix = f"{prefix}.{k_str}" if prefix else k_str
                walk(v, new_prefix)

        elif isinstance(x, list):
            for i, v in enumerate(x):
                new_prefix = f"{prefix}[{i}]" if prefix else f"[{i}]"
                walk(v, new_prefix)

    walk(obj, "")
    return freq_map


def flatten_floats_with_freq_labels(
    obj: Any,
    freq_map: Dict[Tuple[str, int], float],
    prefix: str = ""
) -> Iterable[Tuple[str, float]]:
    """
    Yield (path, float_value) for every float scalar and float list element in obj.
    For list elements, if there is a stable frequency label available in freq_map for
    (list_path, index), label as: "<list_path>[f=<freq>]" instead of "[<index>]".
    """
    if is_float(obj):
        yield (prefix if prefix else "ROOT", obj)
        return

    if isinstance(obj, dict):
        for k, v in obj.items():
            k_str = str(k)
            new_prefix = f"{prefix}.{k_str}" if prefix else k_str
            yield from flatten_floats_with_freq_labels(v, freq_map, new_prefix)
        return

    if isinstance(obj, list):
        # If this list is a pure float-list, it can be labeled per-element using freq_map
        if all(is_float(v) for v in obj):
            list_path = prefix if prefix else "ROOT_LIST"
            for i, v in enumerate(obj):
                fval = freq_map.get((list_path, i))
                if fval is not None:
                    new_prefix = f"{list_path}[f={fval:.6g}]"
                else:
                    new_prefix = f"{list_path}[{i}]"
                yield (new_prefix, v)
            return

        # Otherwise recurse into elements with ordinal paths
        for i, v in enumerate(obj):
            new_prefix = f"{prefix}[{i}]" if prefix else f"[{i}]"
            yield from flatten_floats_with_freq_labels(v, freq_map, new_prefix)
        return

    # ignore everything else


def find_json_files(root_dir: str, filename: str) -> List[str]:
    hits = []
    for dirpath, _, files in os.walk(root_dir):
        if filename in files:
            hits.append(os.path.join(dirpath, filename))
    return sorted(hits)


def load_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def chunked(lst: List[Any], n: int) -> Iterable[List[Any]]:
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def plot_histograms(
    data_by_path: Dict[str, List[float]],
    out_prefix: str,
    max_cols: int,
    max_rows: int,
    bins: int,
    dpi: int,
):
    paths = sorted(data_by_path.keys())
    plots_per_page = max_cols * max_rows

    if not paths:
        print("No floating-point parameters found.")
        return

    page = 1
    for page_paths in chunked(paths, plots_per_page):
        n = len(page_paths)
        cols = min(max_cols, n)
        rows = math.ceil(n / cols)

        fig_w = 4.0 * cols
        fig_h = 3.2 * rows
        fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h))
        if rows == 1 and cols == 1:
            axes = [[axes]]
        elif rows == 1:
            axes = [axes]
        elif cols == 1:
            axes = [[ax] for ax in axes]

        ax_list = [ax for row_axes in axes for ax in row_axes]

        for idx, p in enumerate(page_paths):
            ax = ax_list[idx]
            vals = data_by_path[p]

            ax.hist(vals, bins=bins)
            ax.set_title(p, fontsize=9)
            ax.tick_params(axis="both", labelsize=8)

            if vals:
                vmin = min(vals)
                vmax = max(vals)
                ax.set_xlabel(f"n={len(vals)}  min={vmin:.4g}  max={vmax:.4g}", fontsize=8)

        for j in range(n, len(ax_list)):
            ax_list[j].axis("off")

        fig.suptitle("Parameter Histograms Across Subdirectories", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        out_path = f"{out_prefix}_{page:03d}.png"
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        print(f"Wrote {out_path}")

        page += 1


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Aggregate float parameters from ts.dat.p JSON files in subdirectories and plot histogram grids "
                    "with amplitude/phase sign correction and frequency-labeled list items."
    )
    parser.add_argument("root_dir", help="Root directory containing subdirectories with ts.dat.p")
    parser.add_argument("--filename", default=TARGET_FILENAME, help=f"Target filename (default: {TARGET_FILENAME})")
    parser.add_argument("--out-prefix", default=OUTPUT_PREFIX, help=f"Output prefix (default: {OUTPUT_PREFIX})")
    parser.add_argument("--cols", type=int, default=MAX_COLS, help=f"Max columns per page (default: {MAX_COLS})")
    parser.add_argument("--rows", type=int, default=MAX_ROWS, help=f"Max rows per page (default: {MAX_ROWS})")
    parser.add_argument("--bins", type=int, default=BINS, help=f"Histogram bins (default: {BINS})")
    parser.add_argument("--dpi", type=int, default=DPI, help=f"PNG DPI (default: {DPI})")
    args = parser.parse_args()

    files = find_json_files(args.root_dir, args.filename)
    if not files:
        print(f"No files named '{args.filename}' found under: {args.root_dir}")
        return

    print(f"Found {len(files)} file(s). Parsing + correcting amp/phase + aggregating floats by labeled path...")

    data_by_path: Dict[str, List[float]] = defaultdict(list)

    bad_files = 0
    for fp in files:
        try:
            obj = load_json(fp)
        except Exception as e:
            bad_files += 1
            print(f"[WARN] Failed to parse JSON: {fp} ({e})")
            continue

        # 1) Correct amp/phase sign + wrap in-place
        correct_amp_phase_inplace(obj)

        # 2) Build a map to label list indices with stable frequencies
        freq_map = build_freq_label_map(obj)

        # 3) Flatten floats, using frequency labels where available
        for path, val in flatten_floats_with_freq_labels(obj, freq_map):
            data_by_path[path].append(val)

    print(f"Collected {len(data_by_path)} float-path(s). Bad files: {bad_files}")

    plot_histograms(
        data_by_path=data_by_path,
        out_prefix=args.out_prefix,
        max_cols=args.cols,
        max_rows=args.rows,
        bins=args.bins,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
