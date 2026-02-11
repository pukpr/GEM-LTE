#!/usr/bin/env python3
import sys
import argparse

MISSING = -99999

def parse_value(s: str) -> int:
    # Handles whitespace like "  7068"
    return int(s.strip())

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Subtract mean from 2nd semicolon-delimited field, ignoring -99999."
    )
    ap.add_argument("input", nargs="?", default="-", help="Input file (default: stdin)")
    ap.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    ap.add_argument("--precision", type=int, default=6, help="Decimal places in output (default: 6)")
    args = ap.parse_args()

    # Read all lines (so we can compute mean first)
    if args.input == "-" or args.input == "":
        raw_lines = sys.stdin.read().splitlines()
    else:
        with open(args.input, "r", encoding="utf-8") as f:
            raw_lines = f.read().splitlines()

    rows = []
    valid_vals = []

    for line in raw_lines:
        if not line.strip():
            continue

        parts = line.split(";")
        if len(parts) < 2:
            continue  # skip malformed lines

        field1 = float(parts[0].strip()) - 1.0/24.0
        try:
            v = parse_value(parts[1])
        except ValueError:
            continue  # skip lines with non-numeric 2nd field

        rows.append((field1, v))

        if v != MISSING:
            valid_vals.append(v)

    mean = (sum(valid_vals) / len(valid_vals)) if valid_vals else 0.0

    # Write output: field1 adjusted_value
    out_lines = []
    fmt = "{:." + str(args.precision) + "f}"
    for field1, v in rows:
        if v == MISSING:
            adjusted = 0.0
        else:
            adjusted = float(v) - mean
        out_lines.append(f"{fmt.format(field1)}\t{fmt.format(adjusted)}")

    output_text = "\n".join(out_lines) + ("\n" if out_lines else "")

    if args.output == "-" or args.output == "":
        sys.stdout.write(output_text)
    else:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(output_text)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
