#!/usr/bin/env python3
"""
Compare two time series given as "decimal-year value" pairs (whitespace delimited).
Whitespace may be spaces or tabs (or mixed/multiple spaces) — this script robustly parses any
sequence of whitespace between columns. Comments starting with '#' are ignored.

Features:
- Alignment modes: interp (default), intersection, union
- Compute-only minimal mode: --cc
- Debug messages to stderr: --debug
- Amplitude normalization for each series prior to plotting/computing correlation:
    --norm {none,std,max,zscore}   (default: none)
    --norm-for-corr                (if set, correlation is computed on normalized values)
- Outputs: plot (default) or saved file (--out). In compute-only mode the script prints r (and p if available).

Usage examples:
    # Full behavior (plot + annotation)
    python compare_timeseries_decimal_year_whitespace.py seriesA.txt seriesB.txt

    # Compute-only minimal mode (prints r and exits)
    python compare_timeseries_decimal_year_whitespace.py seriesA.txt seriesB.txt --cc

    # Compute-only with normalization applied for plotting but correlation on original data
    python compare_timeseries_decimal_year_whitespace.py a.txt b.txt --cc --norm max

    # Compute-only with normalization applied and used for correlation
    python compare_timeseries_decimal_year_whitespace.py a.txt b.txt --cc --norm zscore --norm-for-corr --debug
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# optional scipy for p-value
try:
    from scipy.stats import pearsonr  # type: ignore
    _have_scipy = True
except Exception:
    _have_scipy = False


def read_decimal_series(path: Path, debug: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """
    Read a whitespace-delimited two-column file: decimal_year value
    - Accepts spaces, tabs, or any run of whitespace as delimiters.
    - Ignores lines starting with '#'.
    - Accepts files with extra columns; only the first two are used.
    Returns (x: float array of decimal years, y: float array of values)
    """
    if debug:
        print(f"[debug] Reading file '{path}' (whitespace-delimited, comments with '#')", file=sys.stderr)
    try:
        df = pd.read_csv(
            path,
            comment="#",
            delim_whitespace=True,
            header=None,
            usecols=[0, 1],
            names=["x", "y"],
            engine="python",
        )
    except Exception as e:
        raise ValueError(f"Error reading '{path}': {e}")

    if df.empty:
        raise ValueError(f"No data found in '{path}'")

    x = df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    if debug:
        print(f"[debug] Read {x_sorted.size} rows from '{path}'", file=sys.stderr)
    return x_sorted, y_sorted


def decimal_year_to_datetime(decimal_year: float) -> pd.Timestamp:
    year = int(np.floor(decimal_year))
    frac = decimal_year - year
    days = frac * 365.25
    ts = pd.Timestamp(f"{year}-01-01") + pd.to_timedelta(days, unit="D")
    return ts


def align_series(x1: np.ndarray, y1: np.ndarray, x2: np.ndarray, y2: np.ndarray, mode: str, debug: bool = False):
    """
    Align two series according to mode. If debug=True, print steps.
    Returns (x_common, y1_common, y2_common)
    """
    if debug:
        print(f"[debug] Aligning series using mode='{mode}'", file=sys.stderr)

    if mode not in ("interp", "intersection", "union"):
        raise ValueError("mode must be one of 'interp', 'intersection', 'union'")

    idx1 = np.argsort(x1); x1 = x1[idx1]; y1 = y1[idx1]
    idx2 = np.argsort(x2); x2 = x2[idx2]; y2 = y2[idx2]

    if mode == "interp":
        if x2.size < 2:
            raise ValueError("Need at least 2 points in series B to interpolate")
        if debug:
            print(f"[debug] Interpolating series B (size {x2.size}) onto A's {x1.size} x-values", file=sys.stderr)
        y2_interp = np.interp(x1, x2, y2, left=np.nan, right=np.nan)
        valid = ~np.isnan(y2_interp)
        if debug:
            print(f"[debug] After interpolation, {np.count_nonzero(valid)} paired points available", file=sys.stderr)
        return x1[valid], y1[valid], y2_interp[valid]

    if mode == "intersection":
        tol = 1e-8
        common_x = []
        y1_list = []
        y2_list = []
        i = 0
        j = 0
        while i < x1.size and j < x2.size:
            if abs(x1[i] - x2[j]) <= tol:
                common_x.append(x1[i])
                y1_list.append(y1[i])
                y2_list.append(y2[j])
                i += 1
                j += 1
            elif x1[i] < x2[j]:
                i += 1
            else:
                j += 1
        if len(common_x) == 0:
            raise ValueError("No exact-intersection points found between the two series")
        if debug:
            print(f"[debug] Found {len(common_x)} exact-intersection paired points", file=sys.stderr)
        return np.array(common_x), np.array(y1_list), np.array(y2_list)

    # union
    if debug:
        print(f"[debug] Building union of x-values (sizes: A={x1.size}, B={x2.size})", file=sys.stderr)
    x_union = np.union1d(x1, x2)
    if x_union.size == 0:
        raise ValueError("No x-values available for union")
    y1_on_union = np.interp(x_union, x1, y1, left=np.nan, right=np.nan)
    y2_on_union = np.interp(x_union, x2, y2, left=np.nan, right=np.nan)
    valid = ~np.isnan(y1_on_union) & ~np.isnan(y2_on_union)
    if not np.any(valid):
        raise ValueError("No overlapping points after union/interpolation")
    if debug:
        print(f"[debug] After union+interp, {np.count_nonzero(valid)} paired points available", file=sys.stderr)
    return x_union[valid], y1_on_union[valid], y2_on_union[valid]


def normalize_array(y: np.ndarray, method: str, label: str = "", debug: bool = False) -> tuple[np.ndarray, dict]:
    """
    Normalize array y according to method.
    Returns (y_normalized, info_dict) where info_dict contains details (scale, offset, method).
    Supported methods:
      - 'none'   : return y unchanged
      - 'std'    : divide by standard deviation (population, ddof=0)
      - 'max'    : divide by max(abs(y))
      - 'zscore' : (y - mean) / std
    """
    info = {"method": method}
    if method == "none":
        info.update({"scale": None, "offset": None})
        if debug:
            print(f"[debug] No normalization for {label}", file=sys.stderr)
        return y.copy(), info

    if y.size == 0:
        info.update({"scale": None, "offset": None})
        if debug:
            print(f"[debug] Empty array for {label}, returning unchanged", file=sys.stderr)
        return y.copy(), info

    if method == "std":
        scale = np.std(y, ddof=0)
        if scale == 0 or np.isnan(scale):
            if debug:
                print(f"[debug] std is zero or NaN for {label}; skipping scale", file=sys.stderr)
            info.update({"scale": scale, "offset": 0})
            return y.copy(), info
        info.update({"scale": float(scale), "offset": 0.0})
        if debug:
            print(f"[debug] Normalizing {label} by std={scale}", file=sys.stderr)
        return y / scale, info

    if method == "max":
        scale = np.max(np.abs(y))
        if scale == 0 or np.isnan(scale):
            if debug:
                print(f"[debug] max abs is zero or NaN for {label}; skipping scale", file=sys.stderr)
            info.update({"scale": scale, "offset": 0})
            return y.copy(), info
        info.update({"scale": float(scale), "offset": 0.0})
        if debug:
            print(f"[debug] Normalizing {label} by max_abs={scale}", file=sys.stderr)
        return y / scale, info

    if method == "zscore":
        mu = float(np.mean(y))
        sigma = float(np.std(y, ddof=0))
        if sigma == 0 or np.isnan(sigma):
            if debug:
                print(f"[debug] std is zero or NaN for {label}; performing zero-mean but no scale", file=sys.stderr)
            info.update({"scale": sigma, "offset": mu})
            return y - mu, info
        info.update({"scale": sigma, "offset": mu})
        if debug:
            print(f"[debug] Z-scoring {label}: mean={mu}, std={sigma}", file=sys.stderr)
        return (y - mu) / sigma, info

    raise ValueError(f"Unknown normalization method: {method}")


def compute_pearson(y1: np.ndarray, y2: np.ndarray, debug: bool = False) -> tuple[float, float | None]:
    if debug:
        print(f"[debug] Computing Pearson correlation for {y1.size} paired samples", file=sys.stderr)
    if y1.size < 2:
        raise ValueError("Need at least two paired samples to compute correlation")
    if _have_scipy:
        r, p = pearsonr(y1, y2)
        if debug:
            print(f"[debug] scipy.stats.pearsonr returned r={r}, p={p}", file=sys.stderr)
        return float(r), float(p)
    else:
        r = np.corrcoef(y1, y2)[0, 1]
        if debug:
            print(f"[debug] numpy.corrcoef returned r={r} (no p-value available)", file=sys.stderr)
        return float(r), None


def plot_series(x: np.ndarray, y1: np.ndarray, y2: np.ndarray, label1: str, label2: str,
                use_datetime: bool = True, out: Path | None = None, fmt: tuple[str, str] | None = None,
                r: float | None = None, p: float | None = None, n: int | None = None, debug: bool = False) -> None:
    if debug:
        print("[debug] Preparing to create plot", file=sys.stderr)

    if use_datetime:
        datetimes = [decimal_year_to_datetime(xx) for xx in x]
        x_plot = pd.to_datetime(datetimes)
    else:
        x_plot = x

    plt.figure(figsize=(11, 5))
    ax = plt.gca()

    fmt1, fmt2 = ("-", "-") if fmt is None else tuple(fmt.split(",")[:2])

    ax.plot(x_plot, y1, fmt1, label=label1, lw=0.5)
    ax.plot(x_plot, y2, fmt2, label=label2, color='red', lw=0.5)

    if r is not None:
        if p is not None and not math.isnan(p):
            annot = f"Pearson r = {r:.3f}\np = {p:.3g}\nn = {n}"
        else:
            annot = f"Pearson r = {r:.3f}\nn = {n}"
        ax.text(0.02, 0.95, annot, transform=ax.transAxes, va="top", ha="left",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="k"))

    ax.set_xlabel("Date" if use_datetime else "Decimal year")
    ax.set_ylabel("Value")
    ax.set_title(f"Comparison: {label1} vs {label2}")
    ax.grid(alpha=0.25)
    ax.legend()
    plt.tight_layout()
    if out:
        plt.savefig(out, dpi=150)
        if debug:
            print(f"[debug] Saved figure to {out}", file=sys.stderr)
        else:
            print(f"Saved figure to {out}")
    else:
        plt.show()
    plt.close()


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Compare two decimal-year value time series and display Pearson correlation on the plot.")
    p.add_argument("file_a", type=Path, help="First series file (decimal-year value)")
    p.add_argument("file_b", type=Path, help="Second series file (decimal-year value)")
    p.add_argument("--mode", choices=("interp", "intersection", "union"), default="interp",
                   help="How to align the two series before computing correlation (default: interp)")
    p.add_argument("--label-a", default=None, help="Label for series A (default: filename stem)")
    p.add_argument("--label-b", default=None, help="Label for series B (default: filename stem)")
    p.add_argument("--date-x", action="store_true", default=True, help="Use datetime x-axis (convert decimal-year to dates). Enabled by default.")
    p.add_argument("--no-date-x", dest="date_x", action="store_false", help="Do not convert decimal years to datetimes; use numeric x-axis.")
    p.add_argument("--out", type=Path, default=None, help="Save the plot to this file instead of showing it.")
    p.add_argument("--fmt", default=None, help='Two matplotlib formats comma-separated for A and B, e.g. "k-,r--"')
    p.add_argument("--cc", action="store_true", help="Compute correlation coefficient only, print results and exit (no plot).")
    p.add_argument("--debug", action="store_true", help="Print debug messages describing steps to stderr.")
    p.add_argument("--norm", choices=("none", "std", "max", "zscore"), default="none",
                   help="Amplitude normalization applied independently to each series before plotting (default: none).")
    p.add_argument("--norm-for-corr", action="store_true",
                   help="If set, correlation is computed on the normalized series; otherwise correlation uses original values (default).")
    args = p.parse_args(argv)

    debug = args.debug

    # Step 1: Read files
    try:
        if debug:
            print(f"[debug] Starting. file_a={args.file_a}, file_b={args.file_b}, mode={args.mode}, compute-only={args.cc}", file=sys.stderr)
        x1, y1 = read_decimal_series(args.file_a, debug=debug)
        x2, y2 = read_decimal_series(args.file_b, debug=debug)
    except Exception as e:
        print(f"Error reading input files: {e}", file=sys.stderr)
        return 2

    label_a = args.label_a or args.file_a.stem
    label_b = args.label_b or args.file_b.stem

    # Step 2: Align
    try:
        if debug:
            print("[debug] Aligning series now...", file=sys.stderr)
        x_common, y1_common, y2_common = align_series(x1, y1, x2, y2, mode=args.mode, debug=debug)
    except Exception as e:
        print(f"Error aligning series: {e}", file=sys.stderr)
        return 3

    # Keep copies of originals for correlation if needed
    y1_orig = y1_common.copy()
    y2_orig = y2_common.copy()

    # Step 3: Normalize (independently) if requested
    try:
        if args.norm != "none":
            if debug:
                print(f"[debug] Applying normalization method '{args.norm}' to each series", file=sys.stderr)
            y1_norm, info1 = normalize_array(y1_common, args.norm, label=label_a, debug=debug)
            y2_norm, info2 = normalize_array(y2_common, args.norm, label=label_b, debug=debug)
        else:
            if debug:
                print("[debug] No normalization requested", file=sys.stderr)
            y1_norm, info1 = y1_common.copy(), {"method": "none", "scale": None}
            y2_norm, info2 = y2_common.copy(), {"method": "none", "scale": None}
    except Exception as e:
        print(f"Error during normalization: {e}", file=sys.stderr)
        return 6

    if debug:
        print(f"[debug] Normalization info for {label_a}: {info1}", file=sys.stderr)
        print(f"[debug] Normalization info for {label_b}: {info2}", file=sys.stderr)

    # Step 4: Compute Pearson (choose normalized or original based on flag)
    try:
        if args.norm_for_corr:
            if debug:
                print("[debug] Computing correlation on normalized data (as requested by --norm-for-corr)", file=sys.stderr)
            r, p = compute_pearson(y1_norm, y2_norm, debug=debug)
        else:
            if debug:
                print("[debug] Computing correlation on original aligned data (not normalized)", file=sys.stderr)
            r, p = compute_pearson(y1_orig, y2_orig, debug=debug)
    except Exception as e:
        print(f"Error computing Pearson correlation: {e}", file=sys.stderr)
        return 4

    n = y1_common.size

    # If compute-only mode, print minimal output and exit
    if args.cc:
        # Print human-friendly single-line result to stdout, debug details already sent to stderr if requested
        if p is not None:
            print(f"r = {r:.6f}, p = {p:.6g}, n = {n}")
        else:
            print(f"r = {r:.6f}, n = {n} (p-value not available)")
        if debug:
            print("[debug] Exiting due to --cc (compute-only) flag", file=sys.stderr)
        return 0

    # Otherwise, produce the plot (full mode) using normalized series for plotting
    try:
        plot_series(x_common, y1_norm, y2_norm, label_a, label_b,
                    use_datetime=args.date_x, out=args.out, fmt=args.fmt, r=r, p=p, n=n, debug=debug)
    except Exception as e:
        print(f"Error plotting series: {e}", file=sys.stderr)
        return 5

    return 0


if __name__ == "__main__":
    raise SystemExit(main())