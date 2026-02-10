#!/usr/bin/env python3
"""
Compare two time series given as "decimal-year value" pairs (whitespace delimited).
Whitespace may be spaces or tabs (or mixed/multiple spaces) — this script robustly parses any
sequence of whitespace between columns. Comments starting with '#' are ignored.

Usage:
    python compare_timeseries_decimal_year_whitespace.py seriesA.txt seriesB.txt
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.stats import pearsonr  # type: ignore
    _have_scipy = True
except Exception:
    _have_scipy = False


def read_decimal_series(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read a whitespace-delimited two-column file: decimal_year value
    - Accepts spaces, tabs, or any run of whitespace as delimiters.
    - Ignores lines starting with '#'.
    - Accepts files with extra columns; only the first two are used.
    Returns (x: float array of decimal years, y: float array of values)
    """
    # use pandas with delim_whitespace=True which treats any run of whitespace (spaces/tabs) as delimiter
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
    return x[order], y[order]


def decimal_year_to_datetime(decimal_year: float) -> pd.Timestamp:
    year = int(np.floor(decimal_year))
    frac = decimal_year - year
    days = frac * 365.25
    ts = pd.Timestamp(f"{year}-01-01") + pd.to_timedelta(days, unit="D")
    return ts


def align_series(x1: np.ndarray, y1: np.ndarray, x2: np.ndarray, y2: np.ndarray, mode: str):
    if mode not in ("interp", "intersection", "union"):
        raise ValueError("mode must be one of 'interp', 'intersection', 'union'")

    idx1 = np.argsort(x1); x1 = x1[idx1]; y1 = y1[idx1]
    idx2 = np.argsort(x2); x2 = x2[idx2]; y2 = y2[idx2]

    if mode == "interp":
        if x2.size < 2:
            raise ValueError("Need at least 2 points in series B to interpolate")
        y2_interp = np.interp(x1, x2, y2, left=np.nan, right=np.nan)
        valid = ~np.isnan(y2_interp)
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
        return np.array(common_x), np.array(y1_list), np.array(y2_list)

    x_union = np.union1d(x1, x2)
    if x_union.size == 0:
        raise ValueError("No x-values available for union")
    y1_on_union = np.interp(x_union, x1, y1, left=np.nan, right=np.nan)
    y2_on_union = np.interp(x_union, x2, y2, left=np.nan, right=np.nan)
    valid = ~np.isnan(y1_on_union) & ~np.isnan(y2_on_union)
    if not np.any(valid):
        raise ValueError("No overlapping points after union/interpolation")
    return x_union[valid], y1_on_union[valid], y2_on_union[valid]


def compute_pearson(y1: np.ndarray, y2: np.ndarray) -> tuple[float, float | None]:
    if y1.size < 2:
        raise ValueError("Need at least two paired samples to compute correlation")
    if _have_scipy:
        r, p = pearsonr(y1, y2)
        return float(r), float(p)
    else:
        r = np.corrcoef(y1, y2)[0, 1]
        return float(r), None


def plot_series(x: np.ndarray, y1: np.ndarray, y2: np.ndarray, label1: str, label2: str,
                use_datetime: bool = True, out: Path | None = None, fmt: tuple[str, str] | None = None,
                r: float | None = None, p: float | None = None, n: int | None = None) -> None:
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
        if p is not None:
            annot = f"Pearson r = {r:.6f}\np = {p:.6g}\nn = {n}"
        else:
            annot = f"Pearson r = {r:.6f}\nn = {n}"
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
    args = p.parse_args(argv)

    try:
        x1, y1 = read_decimal_series(args.file_a)
        x2, y2 = read_decimal_series(args.file_b)
    except Exception as e:
        print(f"Error reading input files: {e}", file=sys.stderr)
        return 2

    label_a = args.label_a or args.file_a.stem
    label_b = args.label_b or args.file_b.stem

    try:
        x_common, y1_common, y2_common = align_series(x1, y1, x2, y2, mode=args.mode)
    except Exception as e:
        print(f"Error aligning series: {e}", file=sys.stderr)
        return 3

    try:
        r, p = compute_pearson(y1_common, y2_common)
    except Exception as e:
        print(f"Error computing Pearson correlation: {e}", file=sys.stderr)
        return 4

    n = y1_common.size
    plot_series(x_common, y1_common, y2_common, label_a, label_b, use_datetime=args.date_x, out=args.out, fmt=args.fmt, r=r, p=p, n=n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())