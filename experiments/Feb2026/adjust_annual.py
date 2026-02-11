#!/usr/bin/env python3
import sys
import argparse
import math

MISSING = -99999

def parse_float(s: str) -> float:
    return float(s.strip())

def parse_int(s: str) -> int:
    return int(s.strip())

def solve_linear_system(A, b):
    """Gaussian elimination with partial pivoting for small systems."""
    n = len(A)
    M = [row[:] + [b_i] for row, b_i in zip(A, b)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-12:
            raise ValueError("Singular/ill-conditioned regression matrix.")
        M[col], M[pivot] = M[pivot], M[col]

        piv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= piv

        for r in range(n):
            if r == col:
                continue
            factor = M[r][col]
            if factor == 0:
                continue
            for j in range(col, n + 1):
                M[r][j] -= factor * M[col][j]

    return [M[i][n] for i in range(n)]

def build_features(t, include_seasonal: bool, include_half_year: bool):
    """
    Always includes linear trend: [1, t]
    Optionally includes: sin/cos(2πt) and sin/cos(4πt)
    """
    f = [1.0, t]
    if include_seasonal:
        two_pi = 2.0 * math.pi
        f += [math.sin(two_pi * t), math.cos(two_pi * t)]
        if include_half_year:
            f += [math.sin(2.0 * two_pi * t), math.cos(2.0 * two_pi * t)]
    return f

def fit_regression(ts, ys, include_seasonal: bool, include_half_year: bool):
    X = [build_features(t, include_seasonal, include_half_year) for t in ts]
    p = len(X[0])

    XtX = [[0.0] * p for _ in range(p)]
    Xty = [0.0] * p

    for xrow, y in zip(X, ys):
        for i in range(p):
            Xty[i] += xrow[i] * y
            xi = xrow[i]
            for j in range(p):
                XtX[i][j] += xi * xrow[j]

    return solve_linear_system(XtX, Xty)

def predict(t, beta, include_seasonal: bool, include_half_year: bool):
    f = build_features(t, include_seasonal, include_half_year)
    return sum(b * x for b, x in zip(beta, f))

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Parse semicolon-delimited data; always remove linear trend; optionally remove 1y and 0.5y sinusoids; then mean-subtract."
    )
    ap.add_argument("input", nargs="?", default="-", help="Input file (default: stdin)")
    ap.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    ap.add_argument("--precision", type=int, default=6, help="Decimal places in output (default: 6)")
    ap.add_argument("--remove-seasonal", action="store_true",
                    help="Also remove 1-year and 1/2-year sinusoids via multiple linear regression.")
    ap.add_argument("--no-half-year", action="store_true",
                    help="With --remove-seasonal, remove only 1-year (omit 1/2-year terms).")
    args = ap.parse_args()

    # Read lines
    if args.input == "-" or args.input == "":
        raw_lines = sys.stdin.read().splitlines()
    else:
        with open(args.input, "r", encoding="utf-8") as f:
            raw_lines = f.read().splitlines()

    rows = []
    valid_t = []
    valid_y = []

    for line in raw_lines:
        if not line.strip():
            continue
        parts = line.split(";")
        if len(parts) < 2:
            continue

        try:
            t = parse_float(parts[0])
            v = parse_int(parts[1])
        except ValueError:
            continue

        rows.append((t, v))
        if v != MISSING:
            valid_t.append(t)
            valid_y.append(float(v))

    include_half = not args.no_half_year
    include_seasonal = bool(args.remove_seasonal)

    # Need at least as many points as parameters (and preferably more)
    n_params = 2 + (2 if include_seasonal else 0) + (2 if (include_seasonal and include_half) else 0)
    beta = None
    if len(valid_y) >= n_params:
        beta = fit_regression(valid_t, valid_y, include_seasonal, include_half)
    else:
        # Fall back: if not enough points, do no regression (just mean-subtract),
        # but still honor missing->0. (You can change this behavior if you want.)
        beta = None

    # Compute residuals for valid rows, then mean of residuals
    residuals_valid = []
    for t, v in rows:
        if v == MISSING:
            continue
        y = float(v)
        if beta is not None:
            y -= predict(t, beta, include_seasonal, include_half)
        residuals_valid.append(y)

    mean = (sum(residuals_valid) / len(residuals_valid)) if residuals_valid else 0.0

    fmt = "{:." + str(args.precision) + "f}"
    out_lines = []
    for t, v in rows:
        if v == MISSING:
            out = 0.0
        else:
            y = float(v)
            if beta is not None:
                y -= predict(t, beta, include_seasonal, include_half)
            out = y - mean
        out_lines.append(f"{fmt.format(t-1.0/24.0)}\t{fmt.format(out)}")

    output_text = "\n".join(out_lines) + ("\n" if out_lines else "")

    if args.output == "-" or args.output == "":
        sys.stdout.write(output_text)
    else:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(output_text)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
