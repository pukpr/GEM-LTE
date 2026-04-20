#!/usr/bin/env python3
"""
Latent-manifold discovery experiment for GEM-LTE Feb2026 runs.

This script builds a configurable neural experiment around the model-generated
series in `experiments/Feb2026/*/lte_results.csv` using:

- column 1: decimal-year time
- column 2: GEM-LTE model output
- geographic metadata from `experiments/Feb2026/ID.yml`

It infers each subdirectory as one of:

- coastal
- open_ocean
- aggregate

The implementation is intentionally lightweight: it uses an autograd-backed
NumPy neural model with

1. a shared temporal encoder over per-index windows,
2. a geography-aware message-passing layer using a proximity matrix,
3. a low-dimensional latent state,
4. a latent dynamics module, and
5. a geo-aware decoder.

The default training objective uses reconstruction, latent dynamics
consistency, and latent smoothness. Spectral entropy is computed and reported as
diagnostic output; the default config leaves its training weight at 0.0 for
this lightweight implementation.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

try:
    import autograd.numpy as anp
    from autograd import grad
except ImportError as exc:  # pragma: no cover - exercised by CLI execution
    raise SystemExit(
        "This experiment requires the Python package 'autograd'. "
        "Create a virtual environment and install it there, e.g.\n\n"
        "  python3 -m venv .venv\n"
        "  .venv/bin/pip install autograd pyyaml pandas matplotlib\n"
    ) from exc


README_OVERRIDES: dict[str, str] = {
    "brestexcl": "MSL site variant for Brest, France. Coastal tide gauge on the Atlantic coast of Brittany.",
}

COASTAL_SITE_OVERRIDES = {"denison", "brestexcl"}
OPEN_OCEAN_SITE_OVERRIDES = {"155", "183"}

TYPE_ORDER = ("coastal", "open_ocean", "aggregate")


@dataclass
class SiteRecord:
    site_id: str
    name: str
    site_class: str
    basin: str
    latitude: float
    longitude: float
    code: str
    start: float | None
    stop: float | None
    country: str
    series_path: Path
    readme_text: str
    record_months: int
    record_years: float


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Train a lightweight geo-aware neural latent-manifold experiment on Feb2026 GEM-LTE outputs."
    )
    parser.add_argument(
        "--config",
        default=str(script_dir / "nn_latent_manifold_config.yml"),
        help="Path to YAML configuration file.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(script_dir / "nn_latent_manifold_outputs"),
        help="Directory where run artifacts will be written.",
    )
    parser.add_argument("--epochs", type=int, help="Override training.epochs from the config.")
    parser.add_argument("--train-stride", type=int, help="Override model.train_stride from the config.")
    parser.add_argument("--seed", type=int, help="Override training.seed from the config.")
    parser.add_argument(
        "--temporal-encoder-mode",
        choices=["tanh", "sin", "hybrid"],
        help="Override model.temporal_encoder_mode from the config.",
    )
    parser.add_argument(
        "--mean-forcing-mode",
        choices=["hint", "composite_residual"],
        help="Enable mean forcing and choose how it is integrated.",
    )
    parser.add_argument(
        "--use-mean-forcing-hint",
        action="store_true",
        help="Enable experiments/Feb2026/mean_forcing.dat as an auxiliary encoder hint.",
    )
    return parser.parse_args()


def resolve_user_path(path_str: str) -> Path:
    path = Path(path_str).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (Path.cwd() / path).resolve()


def load_config(path: Path, args: argparse.Namespace) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if args.epochs is not None:
        cfg["training"]["epochs"] = int(args.epochs)
    if args.train_stride is not None:
        cfg["model"]["train_stride"] = int(args.train_stride)
    if args.seed is not None:
        cfg["training"]["seed"] = int(args.seed)
    if args.temporal_encoder_mode is not None:
        cfg["model"]["temporal_encoder_mode"] = str(args.temporal_encoder_mode)
    if args.mean_forcing_mode is not None:
        cfg["data"].setdefault("mean_forcing_hint", {})["enabled"] = True
        cfg["data"]["mean_forcing_hint"]["mode"] = str(args.mean_forcing_mode)
    if args.use_mean_forcing_hint:
        cfg["data"].setdefault("mean_forcing_hint", {})["enabled"] = True
        cfg["data"]["mean_forcing_hint"]["mode"] = "hint"
    return cfg


def signed_coord(value: Any, hemi: str) -> float:
    mag = float(value)
    if str(hemi).strip().upper() in {"S", "W"}:
        mag *= -1.0
    return mag


def normalize_basin(country: str, readme_text: str, site_id: str) -> str:
    text = f"{country} {readme_text} {site_id}".lower()
    if any(k in text for k in ("baltic", "gulf of finland", "gulf of bothnia", "oresund")):
        return "baltic"
    if "black sea" in text:
        return "black_sea"
    if any(k in text for k in ("mediterranean", "adriatic")):
        return "mediterranean"
    if any(k in text for k in ("gulf of mexico",)):
        return "gulf_of_mexico"
    if any(k in text for k in ("indian ocean", "indian", "iod")):
        return "indian"
    if any(k in text for k in ("pacific", "japan", "hawaii", "enso", "nino", "pdo", "modoki")):
        return "pacific"
    if any(k in text for k in ("arctic",)):
        return "arctic"
    if any(k in text for k in ("se asia", "south china sea")):
        return "se_asia"
    if any(k in text for k in ("southern ocean",)):
        return "southern"
    if any(k in text for k in ("atlantic", "north sea", "irish sea", "english channel", "gulf of maine", "skagerrak", "kattegat")):
        return "atlantic"
    return country.strip().lower().replace(" ", "_") or "unknown"


def infer_site_class(site_id: str, name: str, code: str, readme_text: str) -> str:
    text = f"{name} {readme_text}".lower()
    code_int = None
    try:
        code_int = int(str(code))
    except Exception:
        code_int = None

    if "composite" in text or (code_int is not None and code_int >= 10000):
        return "aggregate"
    if site_id in OPEN_OCEAN_SITE_OVERRIDES or "open-ocean" in text or "open ocean" in text or "central pacific" in text:
        return "open_ocean"
    if site_id.isdigit() or site_id in COASTAL_SITE_OVERRIDES:
        return "coastal"
    if "coastal tide gauge" in text or "tide gauge" in text:
        return "coastal"
    return "open_ocean"


def load_id_map(path: Path) -> dict[str, dict[str, Any]]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError(f"Unsupported ID.yml structure in {path}")
    return {str(k): dict(v) for k, v in raw.items() if isinstance(v, dict)}


def read_site_readme(site_dir: Path) -> str:
    readme = site_dir / "README.md"
    if readme.exists():
        return readme.read_text(encoding="utf-8", errors="replace").strip()
    return README_OVERRIDES.get(site_dir.name, "")


def decimal_years_to_month_index(values: np.ndarray) -> np.ndarray:
    return np.rint(values * 12.0).astype(int)


def month_index_to_decimal_year(values: np.ndarray) -> np.ndarray:
    return values.astype(float) / 12.0


def load_series_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",", usecols=(0, 1), ndmin=2)
    times = data[:, 0]
    values = data[:, 1]
    months = decimal_years_to_month_index(times)
    if len(np.unique(months)) != len(months):
        grouped: dict[int, list[float]] = {}
        for month, value in zip(months, values):
            grouped.setdefault(int(month), []).append(float(value))
        unique_months = np.array(sorted(grouped), dtype=int)
        mean_values = np.array([np.mean(grouped[m]) for m in unique_months], dtype=float)
        return unique_months, mean_values
    return months, values


def load_mean_forcing_hint(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 2:
        raise ValueError(f"Expected at least 2 columns in mean forcing hint file: {path}")
    months = decimal_years_to_month_index(data[:, 0])
    values = data[:, 1].astype(float)
    return months, values


def longest_consecutive_segment(months: list[int]) -> list[int]:
    if not months:
        return []
    best_start = 0
    best_len = 1
    cur_start = 0
    cur_len = 1
    for i in range(1, len(months)):
        if months[i] == months[i - 1] + 1:
            cur_len += 1
        else:
            if cur_len > best_len:
                best_start = cur_start
                best_len = cur_len
            cur_start = i
            cur_len = 1
    if cur_len > best_len:
        best_start = cur_start
        best_len = cur_len
    return months[best_start: best_start + best_len]


def load_site_records(base_dir: Path, id_map: dict[str, dict[str, Any]], cfg: dict[str, Any]) -> list[SiteRecord]:
    exclude_prefixes = tuple(cfg["data"].get("exclude_prefixes", []))
    exclude_sites = set(cfg["data"].get("exclude_sites", []))
    out: list[SiteRecord] = []
    missing_metadata: list[str] = []

    for site_dir in sorted(base_dir.iterdir()):
        if not site_dir.is_dir():
            continue
        site_id = site_dir.name
        if site_id in exclude_sites or site_id.startswith(exclude_prefixes):
            continue
        csv_path = site_dir / "lte_results.csv"
        if cfg["data"].get("require_lte_results", True) and not csv_path.exists():
            continue

        entry = id_map.get(site_id)
        if entry is None:
            missing_metadata.append(site_id)
            continue

        readme_text = read_site_readme(site_dir)
        months, _values = load_series_csv(csv_path)
        latitude = signed_coord(entry.get("Latitude", 0.0), entry.get("N/S", "N"))
        longitude = signed_coord(entry.get("Longitude", 0.0), entry.get("E/W", "E"))
        basin = normalize_basin(str(entry.get("Country", "")), readme_text, site_id)
        name = str(entry.get("Name", site_id))
        code = str(entry.get("Code", ""))
        out.append(
            SiteRecord(
                site_id=site_id,
                name=name,
                site_class=infer_site_class(site_id, name, code, readme_text),
                basin=basin,
                latitude=latitude,
                longitude=longitude,
                code=code,
                start=float(entry["Start"]) if entry.get("Start") is not None else None,
                stop=float(entry["Stop"]) if entry.get("Stop") is not None else None,
                country=str(entry.get("Country", "")),
                series_path=csv_path,
                readme_text=readme_text,
                record_months=int(len(months)),
                record_years=float(len(months) / 12.0),
            )
        )

    if missing_metadata:
        print("Skipping sites with missing metadata:", ", ".join(sorted(missing_metadata)))

    if not out:
        raise ValueError("No usable site records found after applying exclusions.")
    return out


def build_site_weights(records: list[SiteRecord], power: float) -> np.ndarray:
    lengths = np.array([max(record.record_months, 1) for record in records], dtype=float)
    if power == 0.0:
        weights = np.ones_like(lengths)
    else:
        weights = lengths ** power
    weights = weights / np.mean(weights)
    return weights


def build_aligned_matrix(records: list[SiteRecord]) -> tuple[np.ndarray, np.ndarray]:
    series_by_site: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    start_months: list[int] = []
    stop_months: list[int] = []

    for record in records:
        months, values = load_series_csv(record.series_path)
        series_by_site[record.site_id] = (months, values)
        start_months.append(int(months.min()))
        stop_months.append(int(months.max()))

    shared_start = max(start_months)
    shared_stop = min(stop_months)
    if shared_stop <= shared_start:
        raise ValueError("No common overlapping monthly span across the selected sites.")

    common_months = np.arange(shared_start, shared_stop + 1, dtype=int)
    matrix = np.zeros((len(common_months), len(records)), dtype=float)
    for col, record in enumerate(records):
        months, values = series_by_site[record.site_id]
        matrix[:, col] = np.interp(common_months, months, values)

    return month_index_to_decimal_year(common_months), matrix


def align_hint_series(
    times: np.ndarray,
    hint_cfg: dict[str, Any],
    repo_root: Path,
) -> tuple[np.ndarray | None, np.ndarray | None, dict[str, Any]]:
    if not bool(hint_cfg.get("enabled", False)):
        return None, None, {"enabled": False, "mode": "disabled"}

    hint_path = Path(hint_cfg["path"])
    hint_path = hint_path.resolve() if hint_path.is_absolute() else (repo_root / hint_path).resolve()
    months, values = load_mean_forcing_hint(hint_path)
    target_months = decimal_years_to_month_index(times)
    if target_months[0] < months.min() or target_months[-1] > months.max():
        raise ValueError(
            f"mean_forcing.dat does not fully cover the shared experiment span:\n"
            f"  file span {month_index_to_decimal_year(np.array([months.min(), months.max()]))[0]:.3f}"
            f" to {month_index_to_decimal_year(np.array([months.min(), months.max()]))[1]:.3f}\n"
            f"  required {times[0]:.3f} to {times[-1]:.3f}"
        )
    aligned = np.interp(target_months, months, values)
    mean = float(aligned.mean())
    std = float(aligned.std())
    if std < 1e-8:
        std = 1.0
    normalized = (aligned - mean) / std
    return aligned, normalized, {
        "enabled": True,
        "mode": str(hint_cfg.get("mode", "hint")),
        "path": str(hint_path),
        "mean": mean,
        "std": std,
    }


def build_window_end_indices(n_times: int, window: int, stride: int, require_next: bool) -> np.ndarray:
    max_end = n_times - 2 if require_next else n_times - 1
    return np.arange(window - 1, max_end + 1, stride, dtype=int)


def build_mean_forcing_basis(signal: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            np.ones_like(signal),
            signal,
            signal ** 2,
            signal ** 3,
            np.sin(signal),
            np.cos(signal),
            signal * np.sin(signal),
            signal * np.cos(signal),
        ]
    )


def fit_mean_forcing_baseline(
    signal_raw: np.ndarray,
    signal_norm: np.ndarray,
    target_matrix: np.ndarray,
    fit_indices: np.ndarray,
    ridge: float,
) -> tuple[np.ndarray, np.ndarray]:
    basis_all = build_mean_forcing_basis(signal_norm)
    x_fit = basis_all[fit_indices]
    y_fit = target_matrix[fit_indices]
    gram = x_fit.T @ x_fit + ridge * np.eye(x_fit.shape[1])
    rhs = x_fit.T @ y_fit
    coeffs = np.linalg.solve(gram, rhs)
    baseline = basis_all @ coeffs
    return baseline, coeffs


def standardize_matrix(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means = matrix.mean(axis=0)
    stds = matrix.std(axis=0)
    stds = np.where(stds < 1e-8, 1.0, stds)
    normalized = (matrix - means) / stds
    return normalized, means, stds


def one_hot(index: int, size: int) -> np.ndarray:
    vec = np.zeros(size, dtype=float)
    vec[index] = 1.0
    return vec


def build_geo_features(records: list[SiteRecord]) -> tuple[np.ndarray, list[str]]:
    basins = sorted({record.basin for record in records})
    basin_to_idx = {basin: i for i, basin in enumerate(basins)}
    site_class_to_idx = {name: i for i, name in enumerate(TYPE_ORDER)}

    latitudes = np.array([record.latitude for record in records], dtype=float)
    longitudes = np.array([record.longitude for record in records], dtype=float)
    lat_scale = max(np.abs(latitudes).max(), 1.0)
    lon_scale = max(np.abs(longitudes).max(), 1.0)

    features = []
    for record in records:
        basin_vec = one_hot(basin_to_idx[record.basin], len(basins))
        type_vec = one_hot(site_class_to_idx[record.site_class], len(TYPE_ORDER))
        features.append(
            np.concatenate(
                [
                    np.array([record.latitude / lat_scale, record.longitude / lon_scale], dtype=float),
                    basin_vec,
                    type_vec,
                ]
            )
        )
    return np.vstack(features), basins


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius = 6371.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)
    a = math.sin(d_phi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2.0) ** 2
    return 2.0 * radius * math.atan2(math.sqrt(a), math.sqrt(max(1.0 - a, 0.0)))


def build_proximity_matrix(records: list[SiteRecord], sigma_km: float, basin_block: bool) -> np.ndarray:
    n = len(records)
    prox = np.zeros((n, n), dtype=float)
    for i, left in enumerate(records):
        for j, right in enumerate(records):
            d_km = haversine_km(left.latitude, left.longitude, right.latitude, right.longitude)
            weight = math.exp(-(d_km ** 2) / (2.0 * sigma_km ** 2))
            if basin_block and left.basin != right.basin:
                weight = 0.0
            if i == j:
                weight = 1.0
            prox[i, j] = weight
    row_sums = prox.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0.0, 1.0, row_sums)
    return prox / row_sums


def build_windows(matrix: np.ndarray, times: np.ndarray, window: int, stride: int, require_next: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    xs = []
    ys = []
    ts = []
    max_end = len(times) - 2 if require_next else len(times) - 1
    for end in range(window - 1, max_end + 1, stride):
        xs.append(matrix[end - window + 1: end + 1].T)
        ys.append(matrix[end])
        ts.append(times[end])
    x = np.stack(xs, axis=0)
    y = np.stack(ys, axis=0)
    t = np.array(ts, dtype=float)
    return x, y, t


def build_hint_windows(hint: np.ndarray | None, times: np.ndarray, window: int, stride: int, require_next: bool) -> np.ndarray | None:
    if hint is None:
        return None
    xs = []
    max_end = len(times) - 2 if require_next else len(times) - 1
    for end in range(window - 1, max_end + 1, stride):
        xs.append(hint[end - window + 1: end + 1])
    return np.stack(xs, axis=0)


def split_train_val(
    x: np.ndarray,
    y: np.ndarray,
    t: np.ndarray,
    hint_x: np.ndarray | None,
    fraction: float,
) -> tuple[
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None],
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None],
]:
    n = len(x)
    val_n = max(2, int(round(n * fraction)))
    val_n = min(val_n, n - 2)
    split = n - val_n
    train_hint = None if hint_x is None else hint_x[:split]
    val_hint = None if hint_x is None else hint_x[split:]
    return (x[:split], y[:split], t[:split], train_hint), (x[split:], y[split:], t[split:], val_hint)


def init_params(rng: np.random.Generator, window: int, geo_dim: int, cfg: dict[str, Any]) -> dict[str, anp.ndarray]:
    hidden = cfg["model"]["temporal_hidden_dim"]
    temporal_feature_dim = cfg["model"]["temporal_feature_dim"]
    temporal_encoder_mode = str(cfg["model"].get("temporal_encoder_mode", "tanh")).lower()
    hint_hidden_dim = cfg["model"]["hint_hidden_dim"]
    hint_feature_dim = cfg["model"]["hint_feature_dim"]
    context_dim = cfg["model"]["context_dim"]
    decoder_hidden_dim = cfg["model"]["decoder_hidden_dim"]
    latent_dim = cfg["model"]["latent_dim"]
    temporal_stage2_in_dim = hidden * 2 if temporal_encoder_mode == "hybrid" else hidden
    ext_dim = temporal_feature_dim + geo_dim
    dec_in_dim = latent_dim + geo_dim
    encoder_in_dim = context_dim + hint_feature_dim

    def randn(shape: tuple[int, ...], scale: float = 0.08) -> anp.ndarray:
        return anp.array(rng.normal(0.0, scale, size=shape))

    return {
        "temporal_w1": randn((window, hidden)),
        "temporal_b1": anp.zeros((hidden,)),
        "temporal_w2": randn((temporal_stage2_in_dim, temporal_feature_dim)),
        "temporal_b2": anp.zeros((temporal_feature_dim,)),
        "hint_w1": randn((window, hint_hidden_dim)),
        "hint_b1": anp.zeros((hint_hidden_dim,)),
        "hint_w2": randn((hint_hidden_dim, hint_feature_dim)),
        "hint_b2": anp.zeros((hint_feature_dim,)),
        "context_w_self": randn((ext_dim, context_dim)),
        "context_w_neigh": randn((ext_dim, context_dim)),
        "context_b": anp.zeros((context_dim,)),
        "encoder_w": randn((encoder_in_dim, latent_dim)),
        "encoder_b": anp.zeros((latent_dim,)),
        "dyn_w1": randn((latent_dim, decoder_hidden_dim)),
        "dyn_b1": anp.zeros((decoder_hidden_dim,)),
        "dyn_w2": randn((decoder_hidden_dim, latent_dim)),
        "dyn_b2": anp.zeros((latent_dim,)),
        "decoder_w1": randn((dec_in_dim, decoder_hidden_dim)),
        "decoder_b1": anp.zeros((decoder_hidden_dim,)),
        "decoder_w2": randn((decoder_hidden_dim, 1)),
        "decoder_b2": anp.zeros((1,)),
    }


def encode_windows(
    params: dict[str, anp.ndarray],
    x: anp.ndarray,
    geo: anp.ndarray,
    proximity: anp.ndarray,
    hint_x: anp.ndarray | None = None,
    temporal_encoder_mode: str = "tanh",
) -> anp.ndarray:
    temporal_linear = anp.einsum("nmt,th->nmh", x, params["temporal_w1"]) + params["temporal_b1"]
    if temporal_encoder_mode == "sin":
        temporal_h1 = anp.sin(temporal_linear)
    elif temporal_encoder_mode == "hybrid":
        temporal_h1 = anp.concatenate([anp.tanh(temporal_linear), anp.sin(temporal_linear)], axis=2)
    else:
        temporal_h1 = anp.tanh(temporal_linear)
    temporal_h2 = anp.tanh(anp.einsum("nmh,hf->nmf", temporal_h1, params["temporal_w2"]) + params["temporal_b2"])
    geo_b = anp.broadcast_to(geo[None, :, :], (x.shape[0], geo.shape[0], geo.shape[1]))
    ext = anp.concatenate([temporal_h2, geo_b], axis=2)
    self_term = anp.einsum("nmd,dc->nmc", ext, params["context_w_self"])
    neigh = anp.einsum("ij,njd->nid", proximity, ext)
    neigh_term = anp.einsum("nid,dc->nic", neigh, params["context_w_neigh"])
    context = anp.tanh(self_term + neigh_term + params["context_b"])
    pooled = anp.mean(context, axis=1)
    if hint_x is None:
        hint_features = anp.zeros((x.shape[0], params["hint_b2"].shape[0]))
    else:
        hint_h1 = anp.tanh(anp.dot(hint_x, params["hint_w1"]) + params["hint_b1"])
        hint_features = anp.tanh(anp.dot(hint_h1, params["hint_w2"]) + params["hint_b2"])
    encoder_input = anp.concatenate([pooled, hint_features], axis=1)
    return anp.tanh(anp.dot(encoder_input, params["encoder_w"]) + params["encoder_b"])


def decode_latents(params: dict[str, anp.ndarray], z: anp.ndarray, geo: anp.ndarray) -> anp.ndarray:
    n, latent_dim = z.shape
    m, geo_dim = geo.shape
    z_b = anp.broadcast_to(z[:, None, :], (n, m, latent_dim))
    geo_b = anp.broadcast_to(geo[None, :, :], (n, m, geo_dim))
    dec_in = anp.concatenate([z_b, geo_b], axis=2)
    hidden = anp.tanh(anp.einsum("nmd,dh->nmh", dec_in, params["decoder_w1"]) + params["decoder_b1"])
    return anp.squeeze(anp.einsum("nmh,hq->nmq", hidden, params["decoder_w2"]) + params["decoder_b2"], axis=2)


def latent_dynamics(params: dict[str, anp.ndarray], z: anp.ndarray) -> anp.ndarray:
    hidden = anp.tanh(anp.dot(z, params["dyn_w1"]) + params["dyn_b1"])
    return anp.dot(hidden, params["dyn_w2"]) + params["dyn_b2"]


def loss_components(
    params: dict[str, anp.ndarray],
    x: anp.ndarray,
    y: anp.ndarray,
    geo: anp.ndarray,
    proximity: anp.ndarray,
    hint_x: anp.ndarray | None,
    site_weights: anp.ndarray,
    temporal_encoder_mode: str,
    lambda_dyn: float,
    lambda_smooth: float,
) -> dict[str, anp.ndarray]:
    z = encode_windows(params, x, geo, proximity, hint_x=hint_x, temporal_encoder_mode=temporal_encoder_mode)
    recon = decode_latents(params, z, geo)
    recon_sqerr = (recon - y) ** 2
    recon_loss = anp.sum(recon_sqerr * site_weights[None, :]) / (recon_sqerr.shape[0] * anp.sum(site_weights))

    if len(z) > 1:
        dyn_pred = latent_dynamics(params, z[:-1])
        dyn_loss = anp.mean((dyn_pred - z[1:]) ** 2)
        smooth_loss = anp.mean((z[1:] - z[:-1]) ** 2)
    else:
        dyn_loss = anp.array(0.0)
        smooth_loss = anp.array(0.0)

    total = recon_loss + lambda_dyn * dyn_loss + lambda_smooth * smooth_loss
    return {
        "loss": total,
        "recon_loss": recon_loss,
        "dyn_loss": dyn_loss,
        "smooth_loss": smooth_loss,
    }


def objective(
    params: dict[str, anp.ndarray],
    x: anp.ndarray,
    y: anp.ndarray,
    geo: anp.ndarray,
    proximity: anp.ndarray,
    hint_x: anp.ndarray | None,
    site_weights: anp.ndarray,
    temporal_encoder_mode: str,
    lambda_dyn: float,
    lambda_smooth: float,
) -> anp.ndarray:
    return loss_components(params, x, y, geo, proximity, hint_x, site_weights, temporal_encoder_mode, lambda_dyn, lambda_smooth)["loss"]


def to_float_dict(items: dict[str, Any]) -> dict[str, float]:
    return {key: float(value) for key, value in items.items()}


def adam_step(
    params: dict[str, anp.ndarray],
    grads: dict[str, anp.ndarray],
    m: dict[str, anp.ndarray],
    v: dict[str, anp.ndarray],
    step: int,
    learning_rate: float,
    beta1: float = 0.9,
    beta2: float = 0.999,
    eps: float = 1e-8,
) -> tuple[dict[str, anp.ndarray], dict[str, anp.ndarray], dict[str, anp.ndarray]]:
    new_params: dict[str, anp.ndarray] = {}
    new_m: dict[str, anp.ndarray] = {}
    new_v: dict[str, anp.ndarray] = {}

    for key in params:
        new_m[key] = beta1 * m[key] + (1.0 - beta1) * grads[key]
        new_v[key] = beta2 * v[key] + (1.0 - beta2) * (grads[key] ** 2)
        m_hat = new_m[key] / (1.0 - beta1 ** step)
        v_hat = new_v[key] / (1.0 - beta2 ** step)
        new_params[key] = params[key] - learning_rate * m_hat / (anp.sqrt(v_hat) + eps)
    return new_params, new_m, new_v


def train_model(
    params: dict[str, anp.ndarray],
    train_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None],
    val_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None],
    geo: np.ndarray,
    proximity: np.ndarray,
    site_weights: np.ndarray,
    temporal_encoder_mode: str,
    cfg: dict[str, Any],
) -> tuple[dict[str, anp.ndarray], list[dict[str, float]]]:
    x_train, y_train, _, hint_train = train_data
    x_val, y_val, _, hint_val = val_data
    lambda_dyn = float(cfg["training"]["lambda_dyn"])
    lambda_smooth = float(cfg["training"]["lambda_smooth"])
    loss_grad = grad(objective)
    m = {key: anp.zeros_like(value) for key, value in params.items()}
    v = {key: anp.zeros_like(value) for key, value in params.items()}
    history: list[dict[str, float]] = []
    best_params = copy.deepcopy(params)
    best_val = float("inf")
    patience = 0
    patience_limit = int(cfg["training"]["early_stopping_patience"])
    print_every = int(cfg["training"]["print_every"])
    learning_rate = float(cfg["training"]["learning_rate"])

    for epoch in range(1, int(cfg["training"]["epochs"]) + 1):
        grads = loss_grad(params, x_train, y_train, geo, proximity, hint_train, site_weights, temporal_encoder_mode, lambda_dyn, lambda_smooth)
        params, m, v = adam_step(params, grads, m, v, epoch, learning_rate)
        train_metrics = to_float_dict(loss_components(params, x_train, y_train, geo, proximity, hint_train, site_weights, temporal_encoder_mode, lambda_dyn, lambda_smooth))
        val_metrics = to_float_dict(loss_components(params, x_val, y_val, geo, proximity, hint_val, site_weights, temporal_encoder_mode, lambda_dyn, lambda_smooth))
        row = {"epoch": float(epoch), **{f"train_{k}": v for k, v in train_metrics.items()}, **{f"val_{k}": v for k, v in val_metrics.items()}}
        history.append(row)

        if val_metrics["loss"] < best_val:
            best_val = val_metrics["loss"]
            best_params = copy.deepcopy(params)
            patience = 0
        else:
            patience += 1

        if epoch == 1 or epoch % print_every == 0:
            print(
                f"epoch={epoch:04d} "
                f"train_loss={train_metrics['loss']:.5f} "
                f"val_loss={val_metrics['loss']:.5f} "
                f"recon={val_metrics['recon_loss']:.5f}"
            )

        if patience >= patience_limit:
            print(f"Early stopping at epoch {epoch} (best val loss {best_val:.5f}).")
            break

    return best_params, history


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def top_spectral_peaks(z: np.ndarray, dt_years: float, top_n: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    freqs = np.fft.rfftfreq(z.shape[0], d=dt_years)
    for dim in range(z.shape[1]):
        centered = z[:, dim] - np.mean(z[:, dim])
        spectrum = np.abs(np.fft.rfft(centered)) ** 2
        if len(spectrum) <= 1:
            continue
        valid_freqs = freqs[1:]
        valid_power = spectrum[1:]
        order = np.argsort(valid_power)[::-1][:top_n]
        for rank, idx in enumerate(order, start=1):
            if valid_freqs[idx] <= 0:
                continue
            rows.append(
                {
                    "latent_dim": dim + 1,
                    "rank": rank,
                    "frequency_per_year": float(valid_freqs[idx]),
                    "period_years": float(1.0 / valid_freqs[idx]),
                    "power": float(valid_power[idx]),
                }
            )
    return pd.DataFrame(rows)


def attach_period_matches(peaks: pd.DataFrame, period_library: dict[str, float]) -> pd.DataFrame:
    candidates = [(name, float(period)) for name, period in period_library.items()]
    nearest_names = []
    nearest_periods = []
    rel_diffs = []
    for _, row in peaks.iterrows():
        best_name, best_period = min(candidates, key=lambda item: abs(row["period_years"] - item[1]))
        nearest_names.append(best_name)
        nearest_periods.append(best_period)
        rel_diffs.append(abs(row["period_years"] - best_period) / best_period)
    peaks = peaks.copy()
    peaks["nearest_candidate"] = nearest_names
    peaks["candidate_period_years"] = nearest_periods
    peaks["relative_difference"] = rel_diffs
    return peaks


def spectral_entropy(z: np.ndarray) -> pd.DataFrame:
    rows = []
    for dim in range(z.shape[1]):
        centered = z[:, dim] - np.mean(z[:, dim])
        amps = np.abs(np.fft.rfft(centered))[1:]
        if len(amps) == 0 or np.sum(amps) == 0.0:
            entropy = 0.0
        else:
            probs = amps / np.sum(amps)
            entropy = float(-(probs * np.log(probs + 1e-12)).sum())
        rows.append({"latent_dim": dim + 1, "spectral_entropy": entropy})
    return pd.DataFrame(rows)


def write_metadata(records: list[SiteRecord], output_dir: Path) -> pd.DataFrame:
    df = pd.DataFrame(
        [
            {
                "site_id": record.site_id,
                "name": record.name,
                "site_class": record.site_class,
                "basin": record.basin,
                "latitude": record.latitude,
                "longitude": record.longitude,
                "code": record.code,
                "start_year": record.start,
                "stop_year": record.stop,
                "country": record.country,
                "record_months": record.record_months,
                "record_years": record.record_years,
            }
            for record in records
        ]
    )
    df.to_csv(output_dir / "site_metadata.csv", index=False)
    return df


def plot_latent_spectra(z: np.ndarray, dt_years: float, output_path: Path) -> None:
    fig, axs = plt.subplots(z.shape[1], 1, figsize=(10, max(3, 2.4 * z.shape[1])), sharex=True)
    if z.shape[1] == 1:
        axs = [axs]
    freqs = np.fft.rfftfreq(z.shape[0], d=dt_years)
    for dim, ax in enumerate(axs):
        centered = z[:, dim] - np.mean(z[:, dim])
        power = np.abs(np.fft.rfft(centered)) ** 2
        ax.plot(freqs[1:], power[1:], color="steelblue", linewidth=1.0)
        ax.set_ylabel(f"z{dim + 1}")
        ax.grid(True, ls=":", lw=0.4, alpha=0.7)
    axs[-1].set_xlabel("Frequency (cycles/year)")
    fig.suptitle("Latent power spectra", fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_reconstruction_examples(
    times: np.ndarray,
    truth: np.ndarray,
    recon: np.ndarray,
    records: list[SiteRecord],
    metrics: pd.DataFrame,
    output_path: Path,
    per_class: int,
) -> None:
    chosen: list[int] = []
    for site_class in TYPE_ORDER:
        class_metrics = metrics[metrics["site_class"] == site_class].sort_values("corr", ascending=False)
        for site_id in class_metrics.head(per_class)["site_id"].tolist():
            chosen.append(next(i for i, record in enumerate(records) if record.site_id == site_id))
    if not chosen:
        return

    fig, axs = plt.subplots(len(chosen), 1, figsize=(12, max(3.2, 2.7 * len(chosen))), sharex=True)
    if len(chosen) == 1:
        axs = [axs]
    for ax, idx in zip(axs, chosen):
        record = records[idx]
        ax.plot(times, truth[:, idx], color="black", linewidth=0.9, label="Observed model output")
        ax.plot(times, recon[:, idx], color="tab:red", linewidth=0.8, alpha=0.85, label="NN reconstruction")
        ax.set_title(f"{record.site_id} — {record.name} ({record.site_class}, {record.basin})", fontsize=10)
        ax.grid(True, ls=":", lw=0.4, alpha=0.7)
    axs[0].legend(loc="upper left")
    axs[-1].set_xlabel("Decimal year")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def summarize_metrics(truth: np.ndarray, pred: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rmse = np.sqrt(np.mean((truth - pred) ** 2, axis=0))
    mae = np.mean(np.abs(truth - pred), axis=0)
    corr = np.array([pearson_corr(truth[:, i], pred[:, i]) for i in range(truth.shape[1])], dtype=float)
    return rmse, mae, corr


def fractional_corr_contribution(component_corr: np.ndarray, total_corr: np.ndarray) -> np.ndarray:
    component_corr = np.asarray(component_corr, dtype=float)
    total_corr = np.asarray(total_corr, dtype=float)
    out = np.full_like(component_corr, np.nan, dtype=float)
    valid = np.abs(total_corr) > 1.0e-8
    out[valid] = component_corr[valid] / total_corr[valid]
    return out


def compare_latent_to_mean_forcing(
    times: np.ndarray,
    mean_signal_raw: np.ndarray | None,
    mean_signal_norm: np.ndarray | None,
    eval_end_indices: np.ndarray,
    z_eval: np.ndarray,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    if mean_signal_raw is None or mean_signal_norm is None:
        return None, None

    mean_eval_raw = mean_signal_raw[eval_end_indices]
    mean_eval_norm = mean_signal_norm[eval_end_indices]
    rows: list[dict[str, Any]] = []
    aligned = pd.DataFrame(
        {
            "time": times,
            "mean_forcing_raw": mean_eval_raw,
            "mean_forcing_normalized": mean_eval_norm,
        }
    )

    best_dim = 0
    best_abs_corr = -1.0
    best_sign = 1.0
    best_series = z_eval[:, 0]
    for dim in range(z_eval.shape[1]):
        latent = z_eval[:, dim]
        corr = pearson_corr(mean_eval_norm, latent)
        slope, intercept = np.polyfit(latent, mean_eval_norm, 1)
        rows.append(
            {
                "latent_dim": dim + 1,
                "corr_with_mean_forcing": corr,
                "abs_corr_with_mean_forcing": abs(corr),
                "slope_mean_on_latent": float(slope),
                "intercept_mean_on_latent": float(intercept),
            }
        )
        aligned[f"z{dim + 1}"] = latent
        if abs(corr) > best_abs_corr:
            best_abs_corr = abs(corr)
            best_dim = dim + 1
            best_sign = 1.0 if corr >= 0.0 else -1.0
            best_series = latent

    aligned["best_latent_dim"] = np.full(len(aligned), best_dim)
    aligned["best_latent_signed"] = best_sign * best_series
    comparison_df = pd.DataFrame(rows).sort_values("abs_corr_with_mean_forcing", ascending=False)
    return comparison_df, aligned


def plot_latent_mean_forcing_comparison(
    aligned_df: pd.DataFrame | None,
    comparison_df: pd.DataFrame | None,
    output_path: Path,
) -> None:
    if aligned_df is None or comparison_df is None or aligned_df.empty or comparison_df.empty:
        return

    times = aligned_df["time"].to_numpy()
    mean_norm = aligned_df["mean_forcing_normalized"].to_numpy()
    latent_cols = [col for col in aligned_df.columns if col.startswith("z")]
    best_dim = int(comparison_df.iloc[0]["latent_dim"])
    best_corr = float(comparison_df.iloc[0]["corr_with_mean_forcing"])
    best_latent = aligned_df[f"z{best_dim}"].to_numpy()
    if best_corr < 0.0:
        best_latent = -best_latent

    fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    axs[0].plot(times, mean_norm, color="black", linewidth=1.2, label="mean_forcing.dat")
    for col in latent_cols:
        series = aligned_df[col].to_numpy()
        series = (series - np.mean(series)) / max(np.std(series), 1.0e-8)
        axs[0].plot(times, series, linewidth=0.8, alpha=0.7, label=col)
    axs[0].set_title("Standardized latent dimensions vs mean_forcing.dat")
    axs[0].set_ylabel("z-score")
    axs[0].grid(True, ls=":", lw=0.4, alpha=0.7)
    axs[0].legend(loc="upper right", ncol=2, fontsize=8)

    axs[1].plot(times, mean_norm, color="black", linewidth=1.2, label="mean_forcing.dat")
    best_latent = (best_latent - np.mean(best_latent)) / max(np.std(best_latent), 1.0e-8)
    axs[1].plot(times, best_latent, color="tab:red", linewidth=1.0, label=f"best latent z{best_dim} (signed)")
    axs[1].set_title(f"Best latent comparison (z{best_dim}, corr={abs(best_corr):.4f})")
    axs[1].set_xlabel("Decimal year")
    axs[1].set_ylabel("z-score")
    axs[1].grid(True, ls=":", lw=0.4, alpha=0.7)
    axs[1].legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_ready(v) for v in value]
    return value


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    config_path = resolve_user_path(args.config)
    output_dir = resolve_user_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_config(config_path, args)
    base_dir = (repo_root / cfg["data"]["base_dir"]).resolve()
    id_file = (repo_root / cfg["data"]["id_file"]).resolve()
    rng = np.random.default_rng(int(cfg["training"]["seed"]))

    id_map = load_id_map(id_file)
    records = load_site_records(base_dir, id_map, cfg)
    write_metadata(records, output_dir)

    times, raw_matrix = build_aligned_matrix(records)
    mean_cfg = cfg["data"].get("mean_forcing_hint", {})
    mean_signal_raw, mean_signal_norm, hint_info = align_hint_series(times, mean_cfg, repo_root)
    geo_features, basin_labels = build_geo_features(records)
    proximity = build_proximity_matrix(records, float(cfg["data"]["sigma_km"]), bool(cfg["data"].get("basin_block", True)))
    site_weights = build_site_weights(records, float(cfg["training"].get("record_length_weight_power", 0.0)))
    temporal_encoder_mode = str(cfg["model"].get("temporal_encoder_mode", "tanh")).lower()
    mean_forcing_mode = str(mean_cfg.get("mode", "hint")).lower() if hint_info.get("enabled") else "disabled"

    window = int(cfg["model"]["window_months"])
    train_stride = int(cfg["model"]["train_stride"])
    inference_stride = int(cfg["model"]["inference_stride"])
    train_end_indices = build_window_end_indices(len(times), window, train_stride, require_next=True)
    eval_end_indices = build_window_end_indices(len(times), window, inference_stride, require_next=False)

    mean_baseline_full = np.zeros_like(raw_matrix)
    baseline_coeffs = None
    model_target_matrix = raw_matrix
    if hint_info.get("enabled") and mean_forcing_mode == "composite_residual":
        ridge = float(mean_cfg.get("composite_ridge", 1.0e-4))
        mean_baseline_full, baseline_coeffs = fit_mean_forcing_baseline(
            mean_signal_raw,
            mean_signal_norm,
            raw_matrix,
            train_end_indices,
            ridge=ridge,
        )
        model_target_matrix = raw_matrix - mean_baseline_full

    norm_matrix, means, stds = standardize_matrix(model_target_matrix)

    x_train_all, y_train_all, t_train_all = build_windows(norm_matrix, times, window, train_stride, require_next=True)
    hint_train_all = None
    if hint_info.get("enabled") and mean_forcing_mode == "hint":
        hint_train_all = build_hint_windows(mean_signal_norm, times, window, train_stride, require_next=True)
    train_data, val_data = split_train_val(
        x_train_all,
        y_train_all,
        t_train_all,
        hint_train_all,
        float(cfg["training"]["validation_fraction"]),
    )

    params = init_params(rng, window, geo_features.shape[1], cfg)
    best_params, history = train_model(params, train_data, val_data, geo_features, proximity, site_weights, temporal_encoder_mode, cfg)

    x_eval, y_eval, t_eval = build_windows(norm_matrix, times, window, inference_stride, require_next=False)
    hint_eval = None
    if hint_info.get("enabled") and mean_forcing_mode == "hint":
        hint_eval = build_hint_windows(mean_signal_norm, times, window, inference_stride, require_next=False)
    z_eval = np.array(encode_windows(best_params, x_eval, geo_features, proximity, hint_x=hint_eval, temporal_encoder_mode=temporal_encoder_mode))
    recon_eval_norm = np.array(decode_latents(best_params, z_eval, geo_features))
    residual_pred_eval = recon_eval_norm * stds[None, :] + means[None, :]
    residual_truth_eval = y_eval * stds[None, :] + means[None, :]

    raw_truth_eval = raw_matrix[eval_end_indices]
    mean_only_eval = mean_baseline_full[eval_end_indices]
    if hint_info.get("enabled") and mean_forcing_mode == "composite_residual":
        recon_eval = residual_pred_eval + mean_only_eval
        truth_eval = raw_truth_eval
    else:
        recon_eval = residual_pred_eval
        truth_eval = raw_truth_eval if mean_forcing_mode == "composite_residual" else residual_truth_eval

    metrics_rows = []
    final_rmse, final_mae, final_corr = summarize_metrics(truth_eval, recon_eval)
    baseline_rmse, baseline_mae, baseline_corr = summarize_metrics(raw_truth_eval, mean_only_eval)
    residual_rmse, residual_mae, residual_corr = summarize_metrics(residual_truth_eval, residual_pred_eval)
    residual_corr_increment = final_corr - baseline_corr
    baseline_corr_fraction = fractional_corr_contribution(baseline_corr, final_corr)
    residual_corr_fraction = fractional_corr_contribution(residual_corr_increment, final_corr)
    for idx, record in enumerate(records):
        metrics_rows.append(
            {
                "site_id": record.site_id,
                "name": record.name,
                "site_class": record.site_class,
                "basin": record.basin,
                "rmse": float(final_rmse[idx]),
                "mae": float(final_mae[idx]),
                "corr": float(final_corr[idx]),
                "baseline_rmse": float(baseline_rmse[idx]),
                "baseline_mae": float(baseline_mae[idx]),
                "baseline_corr": float(baseline_corr[idx]),
                "baseline_corr_fraction_of_final": float(baseline_corr_fraction[idx]),
                "residual_rmse": float(residual_rmse[idx]),
                "residual_mae": float(residual_mae[idx]),
                "residual_corr": float(residual_corr[idx]),
                "residual_corr_increment": float(residual_corr_increment[idx]),
                "residual_corr_fraction_of_final": float(residual_corr_fraction[idx]),
                "record_months": int(record.record_months),
                "record_years": float(record.record_years),
                "training_weight": float(site_weights[idx]),
            }
        )
    metrics_df = pd.DataFrame(metrics_rows).sort_values(["site_class", "corr"], ascending=[True, False])
    metrics_df.to_csv(output_dir / "reconstruction_metrics.csv", index=False)

    latent_df = pd.DataFrame({"time": t_eval})
    for i in range(z_eval.shape[1]):
        latent_df[f"z{i + 1}"] = z_eval[:, i]
    latent_df.to_csv(output_dir / "latent_states.csv", index=False)

    latent_mean_df, latent_mean_aligned_df = compare_latent_to_mean_forcing(
        t_eval,
        mean_signal_raw,
        mean_signal_norm,
        eval_end_indices,
        z_eval,
    )
    if latent_mean_df is not None and latent_mean_aligned_df is not None:
        latent_mean_df.to_csv(output_dir / "latent_mean_forcing_comparison.csv", index=False)
        latent_mean_aligned_df.to_csv(output_dir / "latent_mean_forcing_aligned.csv", index=False)

    aligned_df = pd.DataFrame({"time": times})
    for i, record in enumerate(records):
        aligned_df[record.site_id] = raw_matrix[:, i]
    aligned_df.to_csv(output_dir / "aligned_model_matrix.csv", index=False)

    if mean_signal_norm is not None:
        pd.DataFrame({"time": times, "mean_forcing_raw": mean_signal_raw, "mean_forcing_normalized": mean_signal_norm}).to_csv(
            output_dir / "aligned_mean_forcing_hint.csv",
            index=False,
        )
    if hint_info.get("enabled") and mean_forcing_mode == "composite_residual":
        baseline_df = pd.DataFrame({"time": times})
        for i, record in enumerate(records):
            baseline_df[record.site_id] = mean_baseline_full[:, i]
        baseline_df.to_csv(output_dir / "mean_forcing_baseline_matrix.csv", index=False)

    history_df = pd.DataFrame(history)
    history_df.to_csv(output_dir / "training_history.csv", index=False)

    peaks = top_spectral_peaks(z_eval, dt_years=1.0 / 12.0, top_n=int(cfg["report"]["top_spectral_peaks"]))
    peaks = attach_period_matches(peaks, cfg["report"]["period_library_years"])
    peaks.to_csv(output_dir / "latent_spectral_peaks.csv", index=False)
    spectral_entropy(z_eval).to_csv(output_dir / "latent_spectral_entropy.csv", index=False)

    plot_latent_spectra(z_eval, dt_years=1.0 / 12.0, output_path=output_dir / "latent_spectra.png")
    plot_reconstruction_examples(
        t_eval,
        truth_eval,
        recon_eval,
        records,
        metrics_df,
        output_dir / "reconstruction_examples.png",
        per_class=int(cfg["report"]["example_sites_per_class"]),
    )
    plot_latent_mean_forcing_comparison(
        latent_mean_aligned_df,
        latent_mean_df,
        output_dir / "latent_mean_forcing_comparison.png",
    )

    summary = {
        "config_path": str(config_path),
        "output_dir": str(output_dir),
        "site_count": len(records),
        "site_classes": dict(Counter(record.site_class for record in records)),
        "basins": basin_labels,
        "excluded_prefixes": cfg["data"].get("exclude_prefixes", []),
        "excluded_sites": cfg["data"].get("exclude_sites", []),
        "mean_forcing_hint": hint_info,
        "mean_forcing_mode": mean_forcing_mode,
        "record_length_weight_power": float(cfg["training"].get("record_length_weight_power", 0.0)),
        "temporal_encoder_mode": temporal_encoder_mode,
        "shared_span": {"start_year": float(times[0]), "stop_year": float(times[-1]), "months": int(len(times))},
        "window_months": window,
        "train_windows": int(len(train_data[0])),
        "validation_windows": int(len(val_data[0])),
        "inference_windows": int(len(x_eval)),
        "best_validation_loss": float(history_df["val_loss"].min()),
        "mean_reconstruction_corr": float(metrics_df["corr"].mean()),
        "median_reconstruction_corr": float(metrics_df["corr"].median()),
        "mean_baseline_corr": float(metrics_df["baseline_corr"].mean()),
        "median_baseline_corr": float(metrics_df["baseline_corr"].median()),
        "mean_residual_corr": float(metrics_df["residual_corr"].mean()),
        "median_residual_corr": float(metrics_df["residual_corr"].median()),
        "period_library_years": cfg["report"]["period_library_years"],
    }
    if hint_info.get("enabled") and mean_forcing_mode == "composite_residual":
        mean_final_corr = float(metrics_df["corr"].mean())
        mean_baseline_corr = float(metrics_df["baseline_corr"].mean())
        mean_residual_increment = mean_final_corr - mean_baseline_corr
        summary["corr_contribution_breakdown"] = {
            "combined_mean_corr": mean_final_corr,
            "mean_forcing_mean_corr": mean_baseline_corr,
            "residual_increment_mean_corr": mean_residual_increment,
            "mean_forcing_fraction_of_combined_corr": (
                float(mean_baseline_corr / mean_final_corr) if abs(mean_final_corr) > 1.0e-8 else None
            ),
            "residual_fraction_of_combined_corr": (
                float(mean_residual_increment / mean_final_corr) if abs(mean_final_corr) > 1.0e-8 else None
            ),
        }
    if latent_mean_df is not None and not latent_mean_df.empty:
        best_latent_row = latent_mean_df.iloc[0]
        summary["best_latent_mean_forcing_match"] = {
            "latent_dim": int(best_latent_row["latent_dim"]),
            "corr_with_mean_forcing": float(best_latent_row["corr_with_mean_forcing"]),
            "abs_corr_with_mean_forcing": float(best_latent_row["abs_corr_with_mean_forcing"]),
        }
    (output_dir / "run_summary.json").write_text(json.dumps(json_ready(summary), indent=2), encoding="utf-8")

    print(f"Wrote artifacts to {output_dir}")
    print(f"Sites used: {len(records)}")
    print(f"Shared span: {times[0]:.3f} -> {times[-1]:.3f} ({len(times)} months)")
    print(f"Mean reconstruction correlation: {metrics_df['corr'].mean():.4f}")
    if hint_info.get("enabled") and mean_forcing_mode == "composite_residual":
        breakdown = summary["corr_contribution_breakdown"]
        print(
            "Mean CC breakdown: "
            f"mean_forcing.dat={breakdown['mean_forcing_mean_corr']:.4f} "
            f"({breakdown['mean_forcing_fraction_of_combined_corr']:.1%} of combined), "
            f"residual_NN_increment={breakdown['residual_increment_mean_corr']:.4f} "
            f"({breakdown['residual_fraction_of_combined_corr']:.1%} of combined)"
        )


if __name__ == "__main__":
    main()
