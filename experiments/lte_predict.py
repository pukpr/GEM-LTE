#!/usr/bin/env python3
"""
Deterministic predictor companion to annual_forcing_feature_model.py.

Given a target index subdirectory (or site ID), this script loads the saved
annual forcing-feature model outputs, transfers the learned forcing-basis
coefficients and annual correction parameters to the target site, and generates
the predicted monthly time series over the model's shared span.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

import nn_latent_manifold_experiment as shared


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Predict a site time series from the saved annual forcing-feature model."
    )
    parser.add_argument(
        "index_subdir",
        help="Index subdirectory path or site ID to predict.",
    )
    parser.add_argument(
        "--model-dir",
        default=str(script_dir / "annual_forcing_feature_outputs"),
        help="Directory containing outputs from annual_forcing_feature_model.py.",
    )
    parser.add_argument(
        "--output-csv",
        help="Optional path for the predicted time-series CSV. Defaults to <index_subdir>/lte_predict.csv.",
    )
    parser.add_argument(
        "--summary-json",
        help="Optional path for the prediction summary JSON. Defaults to <index_subdir>/lte_predict_summary.json.",
    )
    parser.add_argument(
        "--plot-path",
        help="Optional path for the prediction plot. Defaults to <index_subdir>/lte_predict.png.",
    )
    parser.add_argument(
        "--show-plot",
        action="store_true",
        help="Display the comparison plot after saving it.",
    )
    parser.add_argument(
        "--target-series",
        choices=["model", "actual"],
        default="model",
        help="Which lte_results.csv series to compare against: column 2 model output or column 3 actual data.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def load_csv_with_site_id(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"site_id": str})


def resolve_target_dir(index_arg: str, base_dir: Path) -> Path:
    raw = Path(index_arg).expanduser()
    if raw.is_absolute():
        return raw.resolve()
    cwd_candidate = (Path.cwd() / raw).resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    base_candidate = (base_dir / raw.name).resolve()
    if base_candidate.exists():
        return base_candidate
    return base_candidate


def resolve_output_paths(args: argparse.Namespace, target_dir: Path) -> tuple[Path, Path, Path]:
    output_csv = shared.resolve_user_path(args.output_csv) if args.output_csv else (target_dir / "lte_predict.csv").resolve()
    summary_json = shared.resolve_user_path(args.summary_json) if args.summary_json else (target_dir / "lte_predict_summary.json").resolve()
    plot_path = shared.resolve_user_path(args.plot_path) if args.plot_path else (target_dir / "lte_predict.png").resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    return output_csv, summary_json, plot_path


def infer_target_record(site_id: str, target_dir: Path, id_map: dict[str, dict[str, Any]]) -> dict[str, Any]:
    entry = id_map.get(site_id)
    if entry is None:
        raise ValueError(f"No metadata found in ID.yml for site '{site_id}'.")
    readme_text = shared.read_site_readme(target_dir)
    name = str(entry.get("Name", site_id))
    code = str(entry.get("Code", ""))
    country = str(entry.get("Country", ""))
    return {
        "site_id": site_id,
        "name": name,
        "site_class": shared.infer_site_class(site_id, name, code, readme_text),
        "basin": shared.normalize_basin(country, readme_text, site_id),
        "latitude": shared.signed_coord(entry.get("Latitude", 0.0), entry.get("N/S", "N")),
        "longitude": shared.signed_coord(entry.get("Longitude", 0.0), entry.get("E/W", "E")),
        "code": code,
        "country": country,
    }


def build_weight_vector(target: dict[str, Any], trained_meta: pd.DataFrame, sigma_km: float, basin_block: bool) -> np.ndarray:
    weights = np.zeros(len(trained_meta), dtype=float)
    for idx, row in trained_meta.iterrows():
        if basin_block and str(row["basin"]) != target["basin"]:
            continue
        d_km = shared.haversine_km(
            float(target["latitude"]),
            float(target["longitude"]),
            float(row["latitude"]),
            float(row["longitude"]),
        )
        weights[idx] = np.exp(-(d_km ** 2) / (2.0 * sigma_km ** 2))

    if np.sum(weights) <= 1.0e-12:
        for idx, row in trained_meta.iterrows():
            d_km = shared.haversine_km(
                float(target["latitude"]),
                float(target["longitude"]),
                float(row["latitude"]),
                float(row["longitude"]),
            )
            weights[idx] = 1.0 / max(d_km, 1.0)

    total = float(np.sum(weights))
    if total <= 1.0e-12:
        raise ValueError("Could not form nonzero transfer weights for target site.")
    return weights / total


def weighted_parameter_row(frame: pd.DataFrame, weights: np.ndarray, value_columns: list[str]) -> dict[str, float]:
    values = frame[value_columns].to_numpy(dtype=float)
    averaged = weights @ values
    return {column: float(averaged[idx]) for idx, column in enumerate(value_columns)}


def build_basis_coeff_lookup(
    model_dir: Path,
    basis_names: list[str],
) -> pd.DataFrame:
    coeff_path = model_dir / "mean_forcing_basis_coefficients.csv"
    if coeff_path.exists():
        return load_csv_with_site_id(coeff_path).set_index("site_id", drop=False)

    detail_path = model_dir / "mean_forcing_basis_determination.csv"
    if not detail_path.exists():
        raise FileNotFoundError(f"Missing basis coefficient artifacts in {model_dir}")
    detail_df = load_csv_with_site_id(detail_path)
    pivot = detail_df.pivot(index="site_id", columns="basis_name", values="basis_coefficient").reset_index()
    meta_cols = detail_df.groupby("site_id")[["name", "site_class", "basin"]].first().reset_index()
    merged = meta_cols.merge(pivot, on="site_id", how="left")
    for basis_name in basis_names:
        if basis_name not in merged.columns:
            merged[basis_name] = 0.0
    return merged.set_index("site_id", drop=False)


def build_annual_parameter_lookup(model_dir: Path) -> pd.DataFrame:
    annual_df = load_csv_with_site_id(model_dir / "annual_site_loadings.csv")
    required = {"annual_initial_level", "annual_increment_mean", "annual_increment_std"}
    if required.issubset(set(annual_df.columns)):
        return annual_df.set_index("site_id", drop=False)

    level_path = model_dir / "annual_level_matrix.csv"
    increment_path = model_dir / "annual_increment_matrix.csv"
    if not level_path.exists() or not increment_path.exists():
        raise FileNotFoundError(
            "annual_site_loadings.csv is missing annual parameter columns and the fallback annual matrices are unavailable."
        )

    level_df = pd.read_csv(level_path)
    increment_df = pd.read_csv(increment_path)
    site_cols = [col for col in annual_df.columns if col.startswith("loading_z")]
    rebuilt = annual_df.copy()
    for idx, row in rebuilt.iterrows():
        site_id = str(row["site_id"])
        if site_id not in level_df.columns or site_id not in increment_df.columns:
            raise KeyError(f"Missing fallback annual series for site '{site_id}' in model outputs.")
        rebuilt.loc[idx, "annual_initial_level"] = float(level_df[site_id].iloc[0])
        rebuilt.loc[idx, "annual_increment_mean"] = float(increment_df[site_id].mean())
        std = float(increment_df[site_id].std(ddof=0))
        rebuilt.loc[idx, "annual_increment_std"] = 1.0 if std < 1.0e-8 else std

    ordered_cols = ["site_id", "name", "site_class", "basin", "annual_initial_level", "annual_increment_mean", "annual_increment_std"]
    ordered_cols.extend(site_cols)
    remaining = [col for col in rebuilt.columns if col not in ordered_cols]
    rebuilt = rebuilt[[*ordered_cols, *remaining]]
    return rebuilt.set_index("site_id", drop=False)


def sort_latent_columns(columns: list[str], prefix: str) -> list[str]:
    return sorted([col for col in columns if col.startswith(prefix)], key=lambda col: int(col[len(prefix):]))


def load_target_series(path: Path, series_mode: str) -> tuple[np.ndarray, np.ndarray]:
    usecol = 1 if series_mode == "model" else 2
    data = np.loadtxt(path, delimiter=",", usecols=(0, usecol), ndmin=2)
    times = data[:, 0]
    values = data[:, 1]
    months = shared.decimal_years_to_month_index(times)
    if len(np.unique(months)) != len(months):
        grouped: dict[int, list[float]] = {}
        for month, value in zip(months, values):
            grouped.setdefault(int(month), []).append(float(value))
        unique_months = np.array(sorted(grouped), dtype=int)
        mean_values = np.array([np.mean(grouped[m]) for m in unique_months], dtype=float)
        return unique_months, mean_values
    return months, values


def predict_series(
    times: np.ndarray,
    mean_signal_norm: np.ndarray,
    annual_latent_df: pd.DataFrame,
    basis_names: list[str],
    harmonic_factors: list[float],
    basis_params: dict[str, float],
    annual_params: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    basis_matrix = shared.build_mean_forcing_basis(mean_signal_norm, harmonic_factors=harmonic_factors)
    coeff_vector = np.array([basis_params.get(name, 0.0) for name in basis_names], dtype=float)
    mean_component = basis_matrix @ coeff_vector

    dz_cols = sort_latent_columns(list(annual_latent_df.columns), "dz")
    loading_cols = [f"loading_z{idx + 1}" for idx in range(len(dz_cols))]
    latent_increments = annual_latent_df[dz_cols].to_numpy(dtype=float)
    loading_vector = np.array([annual_params.get(col, 0.0) for col in loading_cols], dtype=float)
    annual_increment = float(annual_params["annual_increment_mean"]) + float(annual_params["annual_increment_std"]) * (latent_increments @ loading_vector)

    annual_levels = np.zeros(len(annual_latent_df), dtype=float)
    annual_levels[0] = float(annual_params["annual_initial_level"])
    for idx in range(1, len(annual_levels)):
        annual_levels[idx] = annual_levels[idx - 1] + annual_increment[idx]

    years = annual_latent_df["year"].to_numpy(dtype=int)
    year_lookup = {int(year): float(level) for year, level in zip(years, annual_levels)}
    months = shared.decimal_years_to_month_index(times)
    annual_component = np.array([year_lookup[int(month // 12)] for month in months], dtype=float)
    prediction = mean_component + annual_component
    return prediction, mean_component, annual_component


def compare_to_actual(target_dir: Path, pred_times: np.ndarray, prediction: np.ndarray) -> tuple[pd.DataFrame | None, dict[str, float] | None]:
    return compare_to_target_series(target_dir, pred_times, prediction, series_mode="model")


def compare_to_target_series(
    target_dir: Path,
    pred_times: np.ndarray,
    prediction: np.ndarray,
    series_mode: str,
) -> tuple[pd.DataFrame | None, dict[str, float] | None]:
    actual_path = target_dir / "lte_results.csv"
    if not actual_path.exists():
        return None, None
    actual_months, actual_values = load_target_series(actual_path, series_mode)
    pred_months = shared.decimal_years_to_month_index(pred_times)
    overlap_mask = (pred_months >= actual_months.min()) & (pred_months <= actual_months.max())
    if not np.any(overlap_mask):
        return None, None
    overlap_months = pred_months[overlap_mask]
    actual_interp = np.interp(overlap_months, actual_months, actual_values)
    pred_overlap = prediction[overlap_mask]
    rmse, mae, corr = shared.summarize_metrics(actual_interp[:, None], pred_overlap[:, None])
    comparison_df = pd.DataFrame(
        {
            "time": pred_times[overlap_mask],
            "predicted": pred_overlap,
            "actual": actual_interp,
            "error": pred_overlap - actual_interp,
        }
    )
    return comparison_df, {
        "rmse": float(rmse[0]),
        "mae": float(mae[0]),
        "corr": float(corr[0]),
        "points": int(len(comparison_df)),
        "series_mode": series_mode,
    }


def plot_prediction(
    plot_path: Path,
    target: dict[str, Any],
    times: np.ndarray,
    prediction: np.ndarray,
    mean_component: np.ndarray,
    annual_component: np.ndarray,
    comparison_df: pd.DataFrame | None,
    target_series_mode: str,
    show_plot: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(times, prediction, color="tab:red", linewidth=1.1, label="Predicted")
    ax.plot(times, mean_component, color="black", linewidth=0.8, linestyle="--", alpha=0.8, label="mean_forcing component")
    ax.plot(times, annual_component, color="tab:blue", linewidth=0.7, alpha=0.7, label="annual correction")
    if comparison_df is not None and not comparison_df.empty:
        label = "Target actual data" if target_series_mode == "actual" else "Target lte_results model"
        ax.plot(comparison_df["time"], comparison_df["actual"], color="tab:green", linewidth=0.9, label=label)
    ax.set_title(f"LTE prediction for {target['site_id']} — {target['name']}")
    ax.set_xlabel("Decimal year")
    ax.set_ylabel("Model value")
    ax.grid(True, ls=":", lw=0.4, alpha=0.7)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=160, bbox_inches="tight")
    if show_plot:
        plt.show()
    plt.close(fig)


def main() -> None:
    args = parse_args()
    model_dir = shared.resolve_user_path(args.model_dir)
    run_summary_path = model_dir / "run_summary.json"
    if not run_summary_path.exists():
        raise FileNotFoundError(f"Missing run_summary.json in model directory: {model_dir}")

    summary = json.loads(run_summary_path.read_text(encoding="utf-8"))
    config_path = Path(summary["config_path"]).resolve()
    cfg = load_yaml(config_path)
    repo_root = Path(__file__).resolve().parents[1]
    base_dir = (repo_root / cfg["data"]["base_dir"]).resolve()
    id_file = (repo_root / cfg["data"]["id_file"]).resolve()

    target_dir = resolve_target_dir(args.index_subdir, base_dir)
    site_id = target_dir.name if target_dir.name else Path(args.index_subdir).name
    output_csv, summary_json, plot_path = resolve_output_paths(args, target_dir)

    id_map = shared.load_id_map(id_file)
    target = infer_target_record(site_id, target_dir, id_map)

    trained_meta = load_csv_with_site_id(model_dir / "site_metadata.csv")
    trained_meta = trained_meta.set_index("site_id", drop=False)
    basis_names = list(summary["mean_forcing_hint"]["basis_names"])
    harmonic_factors = [float(value) for value in summary["mean_forcing_hint"].get("extra_harmonic_factors", [])]

    basis_df = build_basis_coeff_lookup(model_dir, basis_names)
    annual_df = build_annual_parameter_lookup(model_dir)
    latent_df = pd.read_csv(model_dir / "annual_latent_states.csv")
    mean_hint_df = pd.read_csv(model_dir / "aligned_mean_forcing_hint.csv")

    if site_id in basis_df.index and site_id in annual_df.index:
        source_mode = "exact"
        basis_params = {name: float(basis_df.loc[site_id, name]) for name in basis_names}
        annual_cols = ["annual_initial_level", "annual_increment_mean", "annual_increment_std"] + sort_latent_columns(list(annual_df.columns), "loading_z")
        annual_params = {name: float(annual_df.loc[site_id, name]) for name in annual_cols}
        donor_weights = [{"site_id": site_id, "weight": 1.0}]
    else:
        sigma_km = float(cfg["data"].get("sigma_km", 2000.0))
        basin_block = bool(cfg["data"].get("basin_block", True))
        aligned_meta = trained_meta.loc[basis_df.index.intersection(annual_df.index)].reset_index(drop=True)
        aligned_basis = basis_df.loc[aligned_meta["site_id"]].reset_index(drop=True)
        aligned_annual = annual_df.loc[aligned_meta["site_id"]].reset_index(drop=True)
        weights = build_weight_vector(target, aligned_meta, sigma_km=sigma_km, basin_block=basin_block)
        source_mode = "geo_weighted"
        basis_params = weighted_parameter_row(aligned_basis, weights, basis_names)
        annual_cols = ["annual_initial_level", "annual_increment_mean", "annual_increment_std"] + sort_latent_columns(list(aligned_annual.columns), "loading_z")
        annual_params = weighted_parameter_row(aligned_annual, weights, annual_cols)
        order = np.argsort(weights)[::-1]
        donor_weights = [
            {"site_id": str(aligned_meta.iloc[idx]["site_id"]), "weight": float(weights[idx])}
            for idx in order[: min(8, len(order))]
        ]

    times = mean_hint_df["time"].to_numpy(dtype=float)
    mean_signal_norm = mean_hint_df["mean_forcing_normalized"].to_numpy(dtype=float)
    prediction, mean_component, annual_component = predict_series(
        times,
        mean_signal_norm,
        latent_df,
        basis_names,
        harmonic_factors,
        basis_params,
        annual_params,
    )

    pred_df = pd.DataFrame(
        {
            "time": times,
            "predicted": prediction,
            "mean_forcing_component": mean_component,
            "annual_correction_component": annual_component,
        }
    )
    pred_df.to_csv(output_csv, index=False)

    comparison_df, comparison_metrics = compare_to_target_series(target_dir, times, prediction, series_mode=args.target_series)
    comparison_path = None
    if comparison_df is not None:
        comparison_path = output_csv.with_name("lte_predict_comparison.csv")
        comparison_df.to_csv(comparison_path, index=False)
    plot_prediction(
        plot_path,
        target,
        times,
        prediction,
        mean_component,
        annual_component,
        comparison_df,
        target_series_mode=args.target_series,
        show_plot=bool(args.show_plot),
    )

    summary_payload = {
        "target_site": target,
        "source_mode": source_mode,
        "model_dir": str(model_dir),
        "config_path": str(config_path),
        "prediction_csv": str(output_csv),
        "plot_path": str(plot_path),
        "comparison_csv": str(comparison_path) if comparison_path is not None else None,
        "target_series_mode": str(args.target_series),
        "donor_weights": donor_weights,
        "shared_span": {
            "start_year": float(times[0]),
            "stop_year": float(times[-1]),
            "months": int(len(times)),
        },
        "mean_forcing_hint": summary["mean_forcing_hint"],
    }
    if comparison_metrics is not None:
        summary_payload["comparison_metrics"] = comparison_metrics
    summary_json.write_text(json.dumps(shared.json_ready(summary_payload), indent=2), encoding="utf-8")

    print(f"Wrote prediction to {output_csv}")
    print(f"Wrote plot to {plot_path}")
    print(f"Prediction mode: {source_mode}")
    if comparison_metrics is not None:
        print(
            f"Overlap comparison: corr={comparison_metrics['corr']:.4f} "
            f"rmse={comparison_metrics['rmse']:.4f} "
            f"points={comparison_metrics['points']}"
        )


if __name__ == "__main__":
    main()
