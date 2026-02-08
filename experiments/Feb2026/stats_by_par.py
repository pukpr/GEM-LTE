#!/usr/bin/env python3
import argparse
import statistics
import math
from dataclasses import dataclass
from typing import List, Dict, Tuple

import matplotlib.pyplot as plt


OUTLIER_Z = 3.0   # "major outlier" threshold: |x-mean| > OUTLIER_Z * stdev (per key group)
MIN_N_FOR_OUTLIERS = 8  # don't flag outliers for tiny samples


@dataclass
class Obs:
    file: str
    key_raw: float
    val: float


@dataclass
class Group:
    rep: str
    keys: List[str]
    vals: List[float]
    obs: List[Obs]

TWO_PI = 2.0 * math.pi

AMP_PHASE_PAIRS = {
    "ann1": "ann2",
    "sem1": "sem2"
}


def wrap_0_2pi(x: float) -> float:
    return x % TWO_PI

def fix_amp_phase_pairs_postparse(groups: List[Group]) -> None:
    """
    For each (amp_key -> phase_key):
      within each file, if amp < 0 at index i,
      then amp := |amp| and phase[i] := (phase[i] + π) mod 2π
    Assumes amp/phase rows are ordered consistently per file.
    """
    rep_to_group = {g.rep: g for g in groups}

    for amp_key, phase_key in AMP_PHASE_PAIRS.items():
        g_amp = rep_to_group.get(amp_key)
        g_ph  = rep_to_group.get(phase_key)
        if g_amp is None or g_ph is None:
            continue

        # Build per-file ordered lists of indices
        amp_idx_by_file = {}
        for i, ob in enumerate(g_amp.obs):
            amp_idx_by_file.setdefault(ob.file, []).append(i)

        ph_idx_by_file = {}
        for i, ob in enumerate(g_ph.obs):
            ph_idx_by_file.setdefault(ob.file, []).append(i)

        # Apply correction index-wise
        for f in amp_idx_by_file.keys() & ph_idx_by_file.keys():
            a_idx = amp_idx_by_file[f]
            p_idx = ph_idx_by_file[f]
            n = min(len(a_idx), len(p_idx))

            for j in range(n):
                ia = a_idx[j]
                ip = p_idx[j]
                if g_amp.vals[ia] < 0:
                    g_amp.vals[ia] = abs(g_amp.vals[ia])
                    g_amp.obs[ia] = Obs(
                        file=f, key_raw=amp_key, val=g_amp.vals[ia]
                    )

                    v_adj = wrap_0_2pi(g_ph.vals[ip] + math.pi)
                    g_ph.vals[ip] = v_adj
                    g_ph.obs[ip] = Obs(
                        file=f, key_raw=phase_key, val=v_adj
                    )

def assign_group(groups: List[Group], key: str) -> int:
    for i, g in enumerate(groups):
        if g.rep == key:
            return i
    groups.append(Group(rep=key, keys=[key], vals=[], obs=[]))
    return len(groups) - 1


def iter_input_files(list_file: str) -> List[str]:
    files: List[str] = []
    with open(list_file, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            files.append(s.split()[0])
    return files


def parse_data_file(path: str, groups: List[Group]) -> Tuple[int, int]:
    kept = 0
    skipped = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                skipped += 1
                continue
            try:
                k = parts[0]
                if k == "ltep":
                   skipped += 1
                   continue
                if k == "harm":
                   skipped += 1
                   continue
                v = abs(float(parts[1]))
            except ValueError:
                skipped += 1
                continue

            gi = assign_group(groups, k)
            groups[gi].vals.append(v)
            groups[gi].obs.append(Obs(file=path, key_raw=k, val=v))
            kept += 1
    return kept, skipped


def summarize(groups: List[Group]) -> List[Dict[str, object]]:
    rows = []
    for g in groups:
        n = len(g.vals)
        if n == 0:
            continue
        mean = statistics.fmean(g.vals)
        sd = statistics.stdev(g.vals) if n > 1 else 0.0
        sd_pct = (sd / mean * 100.0) if mean != 0 else float("nan")
        rows.append({
            "key_rep": g.rep,
            "n": n,
            "mean": mean,
            "stdev": sd,
            "stdev_pct_mean": sd_pct,
            "min": min(g.vals),
            "max": max(g.vals),
        })
    rows.sort(key=lambda r: r["key_rep"])
    return rows


def find_major_outliers(groups: List[Group]) -> List[Dict[str, object]]:
    outliers: List[Dict[str, object]] = []
    for g in groups:
        n = len(g.vals)
        if n < MIN_N_FOR_OUTLIERS:
            continue
        mean = statistics.fmean(g.vals)
        sd = statistics.stdev(g.vals) if n > 1 else 0.0
        if sd <= 0.0:
            continue

        thresh = OUTLIER_Z * sd
        for ob in g.obs:
            if abs(ob.val - mean) > thresh:
                z = (ob.val - mean) / sd
                outliers.append({
                    "file": ob.file,
                    "key_rep": g.rep,
                    "key_raw": ob.key_raw,
                    "value": ob.val,
                    "mean": mean,
                    "stdev": sd,
                    "z": z,
                })

    outliers.sort(key=lambda r: abs(r["z"]), reverse=True)
    return outliers


def make_scatter_with_errorbars(groups: List[Group], out_png: str) -> None:
    """
    Scatter all values (Y) for each key-group (X) using categorical X positions,
    labeled by the key representative. Overlay +/- 1σ error bars at the mean.
    """
    xs_all: List[int] = []
    ys_all: List[float] = []
    x_means: List[int] = []
    y_means: List[float] = []
    y_sds: List[float] = []
    x_labels: List[str] = []

    for i, g in enumerate(groups):
        if not g.vals:
            continue

        # All observations at categorical position i
        xs_all.extend([i] * len(g.vals))
        ys_all.extend(g.vals)

        m = statistics.fmean(g.vals)
        sd = statistics.stdev(g.vals) if len(g.vals) > 1 else 0.0

        x_means.append(i)
        y_means.append(m)
        y_sds.append(sd)
        x_labels.append(f"{g.rep}")

    plt.figure(figsize=(max(6, 0.4 * len(x_labels)), 8))

    # Scatter points as dash markers
    plt.plot(xs_all, ys_all, linestyle="None", marker="_", markersize=1)

    # Error bars at group means
    # plt.errorbar(x_means, y_means, yerr=y_sds, fmt="none", capsize=3, ecolor="red")
    plt.errorbar(x_means, y_means, yerr=y_sds, fmt="o", capsize=3, ecolor="red")

    plt.xticks(range(len(x_labels)), x_labels, rotation=45, ha="right")
    plt.xlabel("Key (categorical, group representative)")
    plt.ylabel("Value (field2)")
    plt.axhline(0, linestyle="--", linewidth=0.8)
    # plt.yscale("log")
    plt.title("All values by key with ±1σ error bars")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("list_file", nargs="?", default="list_pars.txt",
                    help="Whitespace-delimited file containing input filenames (default: list_pars.txt)")
    ap.add_argument("--plot", default="scatter_by_par.png",
                    help="Output PNG filename for the scatter+errorbar plot")
    args = ap.parse_args()

    files = iter_input_files(args.list_file)
    if not files:
        raise SystemExit(f"No input files found in {args.list_file}")

    groups: List[Group] = []
    total_kept = 0
    total_skipped = 0

    for fp in files:
        kept, skipped = parse_data_file(fp, groups)
        total_kept += kept
        total_skipped += skipped

    # fix_amp_phase_pairs_postparse(groups)

    stats_rows = summarize(groups)
    outliers = find_major_outliers(groups)

    print(f"# Files: {len(files)}")
    print(f"# Rows kept (exactly 3 fields, numeric key+value): {total_kept}")
    print(f"# Rows skipped: {total_skipped}")
    print()

    print("## Stats by key-group (field1 ~ key, field2 summarized)")
    print(f"{'key(rep)':>14}  {'n':>7}  {'mean(field2)':>14}  {'stdev':>14}  {'stdev%mean':>11}")
    print("-" * 70)
    for r in stats_rows:
        print(f"{r['key_rep']}  {r['n']:7d}  {r['mean']:14.6g}  {r['stdev']:14.6g}  {r['stdev_pct_mean']:11.4g}")

    print("\n## Major outliers")
    print(f"# Rule: n >= {MIN_N_FOR_OUTLIERS}, stdev > 0, and |x-mean| > {OUTLIER_Z}*stdev (per key group)")
    if not outliers:
        print("# None found under current thresholds.")
    else:
        print(f"{'|z|':>8}  {'value':>14}  {'mean':>14}  {'stdev':>14}  {'key(rep)':>14}  {'key(raw)':>14}  file")
        print("-" * 110)
        for o in outliers:
            print(f"{abs(o['z']):8.3f}  {o['value']:14.6g}  {o['mean']:14.6g}  {o['stdev']:14.6g}  "
                  f"{o['key_rep']} {o['file']}")

        files_with_outliers = sorted({o["file"] for o in outliers})
        print("\n# Files featuring major outliers:")
        for f in files_with_outliers:
            print(f"- {f}")

    # Plot
    make_scatter_with_errorbars(groups, args.plot)
    print(f"\n# Wrote plot: {args.plot}")


if __name__ == "__main__":
    main()


