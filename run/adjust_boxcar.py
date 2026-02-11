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
    Optionally adds: sin/cos(2πt) and sin/cos(4πt)
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

def boxcar3_two_pass(values, missing_mask, protect_endpoints=True, exclude_endpoints_from_neighbors=True):
    """
    Two consecutive 3-pt boxcar filters, respecting missing data and endpoints.

    Rules:
    - Endpoints are never changed if protect_endpoints=True.
    - Endpoints do not contribute to neighbors if exclude_endpoints_from_neighbors=True
      (prevents endpoint zeros from dragging down point 1 / n-2).
    - Missing samples are never changed and never used in averages.
    - If only 1 or 2 valid samples exist in the 3-pt window, average over those.
    """
    n = len(values)
    if n < 3:
        return values[:]  # nothing to do

    endpoint_mask = [False] * n
    if protect_endpoints:
        endpoint_mask[0] = True
        endpoint_mask[-1] = True

    def one_pass(x):
        y = x[:]
        for i in range(n):
            if missing_mask[i]:
                y[i] = 0.0
                continue
            if endpoint_mask[i]:
                y[i] = x[i]
                continue

            idxs = [i - 1, i, i + 1]
            acc = 0.0
            cnt = 0

            for j in idxs:
                if j < 0 or j >= n:
                    continue
                if missing_mask[j]:
                    continue
                if exclude_endpoints_from_neighbors and endpoint_mask[j]:
                    continue
                acc += x[j]
                cnt += 1

            # If nothing valid (can happen near endpoints with exclusion), leave as-is
            y[i] = (acc / cnt) if cnt > 0 else x[i]
        return y

    v1 = one_pass(values)
    v2 = one_pass(v1)

    # enforce invariants explicitly
    if protect_endpoints:
        v2[0] = values[0]
        v2[-1] = values[-1]
    for i in range(n):
        if missing_mask[i]:
            v2[i] = 0.0

    return v2

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Parse semicolon-delimited data; always remove linear trend; optionally remove 1y/0.5y sinusoids; mean-subtract; then 2x 3-pt boxcar respecting missing/endpoints."
    )
    ap.add_argument("input", nargs="?", default="-", help="Input file (default: stdin)")
    ap.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    ap.add_argument("--precision", type=int, default=6, help="Decimal places in output (default: 6)")

    ap.add_argument("--remove-seasonal", action="store_true",
                    help="Also remove 1-year and 1/2-year sinusoids via multiple linear regression.")
    ap.add_argument("--no-half-year", action="store_true",
                    help="With --remove-seasonal, remove only 1-year (omit 1/2-year terms).")

    ap.add_argument("--no-boxcar", action="store_true",
                    help="Disable the final two-pass 3-pt boxcar smoothing.")

    args = ap.parse_args()

    # Read input
    if args.input == "-" or args.input == "":
        raw_lines = sys.stdin.read().splitlines()
    else:
        with open(args.input, "r", encoding="utf-8") as f:
            raw_lines = f.read().splitlines()

    years = []
    vals_int = []
    missing_mask = []

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

        years.append(t)
        vals_int.append(v)
        is_missing = (v == MISSING)
        missing_mask.append(is_missing)

        if not is_missing:
            valid_t.append(t)
            valid_y.append(float(v))

    n = len(years)
    if n == 0:
        return 0

    include_half = not args.no_half_year
    include_seasonal = bool(args.remove_seasonal)

    # Fit regression on valid points only
    n_params = 2 + (2 if include_seasonal else 0) + (2 if (include_seasonal and include_half) else 0)
    beta = None
    if len(valid_y) >= n_params:
        beta = fit_regression(valid_t, valid_y, include_seasonal, include_half)

    # Residuals (detrended/deseasoned); missing -> 0 placeholder
    residuals = [0.0] * n
    residuals_valid = []
    for i, (t, v) in enumerate(zip(years, vals_int)):
        if v == MISSING:
            residuals[i] = 0.0
            continue
        y = float(v)
        if beta is not None:
            y -= predict(t, beta, include_seasonal, include_half)
        residuals[i] = y
        residuals_valid.append(y)

    # Subtract mean of valid residuals
    mean = (sum(residuals_valid) / len(residuals_valid)) if residuals_valid else 0.0
    adjusted = [0.0 if missing_mask[i] else (residuals[i] - mean) for i in range(n)]

    # Final smoothing: two consecutive 3-pt boxcars
    if not args.no_boxcar:
        adjusted = boxcar3_two_pass(
            adjusted,
            missing_mask=missing_mask,
            protect_endpoints=True,
            exclude_endpoints_from_neighbors=True
        )

    # Output
    fmt = "{:." + str(args.precision) + "f}"
    out_lines = [f"{fmt.format(years[i]-1.0/24.0)}\t{fmt.format(adjusted[i])}" for i in range(n)]
    output_text = "\n".join(out_lines) + ("\n" if out_lines else "")

    if args.output == "-" or args.output == "":
        sys.stdout.write(output_text)
    else:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(output_text)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
