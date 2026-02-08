#!/usr/bin/env python3
"""
Fill missing monthly steps in a whitespace-delimited time series file, in-place.

Input (whitespace-delimited: spaces/tabs):
  <decimal_year> <value> [optional extra columns...]

Rules:
- Column 1: decimal year (float). Expected monthly step = 1/12 year.
- Column 2: numeric value (float).
- If any consecutive timestamps are separated by more than ~1 month, insert missing
  timestamps at monthly intervals with value 0.0.
- Overwrites the input file safely (temp file + atomic replace).

Notes:
- Uses epsilon (tolerance) with math.isclose to account for float rounding noise.
- Handles gaps of many months efficiently.

Example:
  2025.000000  1.2
  2025.166667  3.4
Will insert:
  2025.083333  0.0
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from dataclasses import dataclass
from typing import List


MONTH_STEP = 1.0 / 12.0


@dataclass(frozen=True)
class Row:
    t: float
    v: float
    extras: List[str]


def is_blank_or_comment(line: str) -> bool:
    s = line.strip()
    return not s or s.startswith("#")


def parse_row(line: str, lineno: int) -> Row:
    parts = line.strip().split()  # whitespace-delimited
    if len(parts) < 2:
        raise ValueError(f"Line {lineno}: expected >=2 columns, got {len(parts)}")
    try:
        t = float(parts[0])
    except ValueError as e:
        raise ValueError(f"Line {lineno}: first column not a float: {parts[0]!r}") from e
    try:
        v = float(parts[1])
    except ValueError as e:
        raise ValueError(f"Line {lineno}: second column not a float: {parts[1]!r}") from e
    return Row(t=t, v=v, extras=parts[2:])


def fmt_time(t: float, decimals: int) -> str:
    return f"{t:.{decimals}f}"


def fmt_value(v: float) -> str:
    # Keep inserted zeros as 0.0, and keep original values compact.
    if v == 0.0:
        return "0.0"
    return f"{v:.6g}"


def is_close(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    return math.isclose(a, b, abs_tol=abs_tol, rel_tol=rel_tol)


def fill_month_gaps(rows: List[Row], abs_tol: float, rel_tol: float) -> List[Row]:
    if not rows:
        return []

    # Sort by time to interpret as a time series.
    rows_sorted = sorted(rows, key=lambda r: r.t)

    out: List[Row] = [rows_sorted[0]]

    for cur in rows_sorted[1:]:
        prev = out[-1]

        # Near-duplicate timestamps: keep cur, don't fill between.
        if is_close(cur.t, prev.t, abs_tol=abs_tol, rel_tol=rel_tol):
            out.append(cur)
            continue

        dt = cur.t - prev.t

        # If dt is ~1 month, just append.
        if is_close(dt, MONTH_STEP, abs_tol=abs_tol, rel_tol=rel_tol):
            out.append(cur)
            continue

        # If dt is larger than one month, insert intermediate months.
        if dt > MONTH_STEP and not is_close(dt, MONTH_STEP, abs_tol=abs_tol, rel_tol=rel_tol):
            # Estimate how many month steps are between prev and cur.
            # Use rounding to handle floating error (e.g., 2.999999 months -> 3).
            approx_steps = int(round(dt / MONTH_STEP))

            # If approx_steps <= 1, nothing to insert.
            if approx_steps <= 1:
                out.append(cur)
                continue

            # Insert month steps: prev + 1*step ... prev + (approx_steps-1)*step,
            # but don't accidentally overshoot/duplicate cur due to rounding.
            for k in range(1, approx_steps):
                t_missing = prev.t + k * MONTH_STEP

                # If this missing time would land on/after cur (within tolerance), stop.
                if (t_missing > cur.t) or is_close(t_missing, cur.t, abs_tol=abs_tol, rel_tol=rel_tol):
                    break

                out.append(Row(t=t_missing, v=0.0, extras=[]))

            out.append(cur)
            continue

        # If dt is less than one month (but not duplicate), keep both; no filling.
        out.append(cur)

    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Fill missing monthly timestamps (1/12 year) with 0.0 values; overwrite file in-place."
    )
    ap.add_argument("filename", help="Input whitespace-delimited file to overwrite.")
    ap.add_argument(
        "--abs-tol",
        type=float,
        default=1e-4 * MONTH_STEP,
        help="Absolute tolerance for float comparisons (default: 1e-4 * (1/12)).",
    )
    ap.add_argument(
        "--rel-tol",
        type=float,
        default=1e-9,
        help="Relative tolerance for float comparisons (default: 1e-9).",
    )
    ap.add_argument(
        "--time-decimals",
        type=int,
        default=6,
        help="Decimals to write for the timestamp column (default: 6).",
    )
    ap.add_argument(
        "--delimiter",
        choices=["tab", "space"],
        default="tab",
        help="Delimiter to use when writing output (default: tab).",
    )
    args = ap.parse_args()

    filename = args.filename
    abs_tol = args.abs_tol
    rel_tol = args.rel_tol
    time_decimals = args.time_decimals
    delim = "\t" if args.delimiter == "tab" else " "

    header_lines: List[str] = []
    data_rows: List[Row] = []

    # Preserve initial blank/comment header lines (common in these files).
    with open(filename, "r", encoding="utf-8") as f:
        saw_data = False
        for lineno, line in enumerate(f, start=1):
            if is_blank_or_comment(line):
                if not saw_data:
                    header_lines.append(line.rstrip("\n"))
                continue
            saw_data = True
            data_rows.append(parse_row(line, lineno))

    filled = fill_month_gaps(data_rows, abs_tol=abs_tol, rel_tol=rel_tol)

    # Write to temp file in same directory, then atomically replace original.
    dirpath = os.path.dirname(os.path.abspath(filename)) or "."
    fd, tmppath = tempfile.mkstemp(prefix=".fill_month_gaps_", dir=dirpath, text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as out:
            for line in header_lines:
                out.write(line + "\n")

            for r in filled:
                cols = [fmt_time(r.t, time_decimals), fmt_value(r.v)]
                if r.extras:
                    cols.extend(r.extras)
                out.write(delim.join(cols) + "\n")

        os.replace(tmppath, filename)
    except Exception:
        try:
            os.unlink(tmppath)
        except OSError:
            pass
        raise

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
