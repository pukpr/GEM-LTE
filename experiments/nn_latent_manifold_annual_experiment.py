#!/usr/bin/env python3
"""
Annualized latent-manifold offshoot for GEM-LTE Feb2026 runs.

This variant keeps the same data sources and most of the same report artifacts as
`nn_latent_manifold_experiment.py`, but replaces the monthly latent NN with an
explicit annual increment factor model:

1. fit the optional nonlinear `mean_forcing.dat` baseline,
2. aggregate the remaining residual to annual means,
3. factor yearly increments in those annual levels with a low-rank latent model,
4. expand the reconstructed annual levels back to monthly resolution.

The goal is to capture sharper year-to-year level changes in an explainable and
regularized way while preserving the experiment outputs used for comparison.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

import nn_latent_manifold_experiment as shared


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Run the annualized latent-manifold offshoot on Feb2026 GEM-LTE outputs."
    )
    parser.add_argument(
        "--config",
        default=str(script_dir / "nn_latent_manifold_annual_config.yml"),
        help="Path to YAML configuration file.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(script_dir / "nn_latent_manifold_annual_outputs"),
        help="Directory where run artifacts will be written.",
    )
    parser.add_argument(
        "--annual-latent-dim",
        type=int,
        help="Override annual_model.latent_dim from the config.",
    )
    return parser.parse_args()


def load_config(path: Path, args: argparse.Namespace) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if args.annual_latent_dim is not None:
        cfg.setdefault("annual_model", {})["latent_dim"] = int(args.annual_latent_dim)
    return cfg


def build_annual_level_matrix(times: np.ndarray, matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    month_indices = shared.decimal_years_to_month_index(times)
    years = month_indices // 12
    annual_years, inverse = np.unique(years, return_inverse=True)
    annual_levels = np.zeros((len(annual_years), matrix.shape[1]), dtype=float)
    months_per_year = np.zeros(len(annual_years), dtype=int)
    for idx in range(len(annual_years)):
        mask = inverse == idx
        months_per_year[idx] = int(mask.sum())
        annual_levels[idx] = matrix[mask].mean(axis=0)
    return annual_years.astype(int), annual_levels, months_per_year, inverse.astype(int)


def fit_annual_increment_factor_model(
    annual_levels: np.ndarray,
    site_weights: np.ndarray,
    latent_dim: int,
) -> dict[str, np.ndarray]:
    if annual_levels.shape[0] < 2:
        raise ValueError("Need at least two annual levels to fit annual increments.")

    annual_increments = np.diff(annual_levels, axis=0)
    increment_mean = annual_increments.mean(axis=0, keepdims=True)
    increment_std = annual_increments.std(axis=0, keepdims=True)
    increment_std = np.where(increment_std < 1.0e-8, 1.0, increment_std)
    increment_norm = (annual_increments - increment_mean) / increment_std

    sqrt_weights = np.sqrt(np.asarray(site_weights, dtype=float))
    weighted = increment_norm * sqrt_weights[None, :]
    u, singular_values, vt = np.linalg.svd(weighted, full_matrices=False)
    rank = max(1, min(int(latent_dim), weighted.shape[0], weighted.shape[1]))

    u_r = u[:, :rank]
    s_r = singular_values[:rank]
    vt_r = vt[:rank, :]
    latent_increments = u_r * s_r[None, :]
    site_loadings = vt_r.T / sqrt_weights[:, None]
    increment_recon_norm = latent_increments @ site_loadings.T
    increment_recon = increment_recon_norm * increment_std + increment_mean

    annual_level_recon = np.zeros_like(annual_levels)
    annual_level_recon[0] = annual_levels[0]
    annual_level_recon[1:] = annual_levels[0] + np.cumsum(increment_recon, axis=0)

    latent_levels = np.zeros((annual_levels.shape[0], rank), dtype=float)
    latent_levels[1:] = np.cumsum(latent_increments, axis=0)

    explained = np.zeros(rank, dtype=float)
    total_energy = float(np.sum(singular_values ** 2))
    if total_energy > 0.0:
        explained = (s_r ** 2) / total_energy

    return {
        "annual_increments": annual_increments,
        "annual_increment_recon": increment_recon,
        "annual_level_recon": annual_level_recon,
        "increment_mean": increment_mean.ravel(),
        "increment_std": increment_std.ravel(),
        "latent_increments": latent_increments,
        "latent_levels": latent_levels,
        "site_loadings": site_loadings,
        "singular_values": s_r,
        "explained_variance_ratio": explained,
    }


def expand_annual_levels(annual_level_matrix: np.ndarray, year_inverse: np.ndarray) -> np.ndarray:
    return annual_level_matrix[year_inverse]


def write_wide_matrix(path: Path, times: np.ndarray, records: list[shared.SiteRecord], matrix: np.ndarray) -> None:
    df = pd.DataFrame({"time": times})
    for idx, record in enumerate(records):
        df[record.site_id] = matrix[:, idx]
    df.to_csv(path, index=False)


def write_annual_latent_states(path: Path, annual_years: np.ndarray, latent_levels: np.ndarray, latent_increments: np.ndarray) -> None:
    df = pd.DataFrame({"year": annual_years})
    for idx in range(latent_levels.shape[1]):
        df[f"z{idx + 1}"] = latent_levels[:, idx]
        increment_col = np.zeros(len(annual_years), dtype=float)
        increment_col[1:] = latent_increments[:, idx]
        df[f"dz{idx + 1}"] = increment_col
    df.to_csv(path, index=False)


def write_site_loadings(path: Path, records: list[shared.SiteRecord], site_loadings: np.ndarray) -> None:
    rows: list[dict[str, Any]] = []
    for site_idx, record in enumerate(records):
        row: dict[str, Any] = {
            "site_id": record.site_id,
            "name": record.name,
            "site_class": record.site_class,
            "basin": record.basin,
        }
        for dim in range(site_loadings.shape[1]):
            row[f"loading_z{dim + 1}"] = float(site_loadings[site_idx, dim])
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def coefficient_of_determination(truth: np.ndarray, pred: np.ndarray) -> np.ndarray:
    truth = np.asarray(truth, dtype=float)
    pred = np.asarray(pred, dtype=float)
    sse = np.sum((truth - pred) ** 2, axis=0)
    centered = truth - np.mean(truth, axis=0, keepdims=True)
    sst = np.sum(centered ** 2, axis=0)
    r2 = np.zeros(truth.shape[1], dtype=float)
    valid = sst > 1.0e-12
    r2[valid] = 1.0 - (sse[valid] / sst[valid])
    return r2


def fit_basis_subset(
    basis_all: np.ndarray,
    target_matrix: np.ndarray,
    keep_indices: list[int],
    ridge: float,
    l1_alpha: float,
) -> tuple[np.ndarray, np.ndarray]:
    x_fit = basis_all[:, keep_indices]
    gram = x_fit.T @ x_fit + ridge * np.eye(x_fit.shape[1])
    rhs = x_fit.T @ target_matrix
    coeffs = np.linalg.solve(gram, rhs)
    pred = x_fit @ coeffs
    return pred, coeffs


def build_mean_forcing_basis_determination(
    basis_all: np.ndarray,
    basis_coeffs: np.ndarray,
    raw_matrix: np.ndarray,
    full_baseline: np.ndarray,
    records: list[shared.SiteRecord],
    basis_names: list[str],
    ridge: float,
    l1_alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_baseline_r2 = coefficient_of_determination(raw_matrix, full_baseline)
    rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for basis_idx, basis_name in enumerate(basis_names):
        single_keep = [0] if basis_idx == 0 else [0, basis_idx]
        single_pred, single_coeffs = fit_basis_subset(basis_all, raw_matrix, single_keep, ridge, l1_alpha)
        single_r2_raw = coefficient_of_determination(raw_matrix, single_pred)
        keep_indices = [idx for idx in range(len(basis_names)) if idx != basis_idx]
        if keep_indices:
            leave_out_pred, _ = fit_basis_subset(basis_all, raw_matrix, keep_indices, ridge, l1_alpha)
            leave_out_r2 = coefficient_of_determination(raw_matrix, leave_out_pred)
        else:
            leave_out_r2 = np.zeros(raw_matrix.shape[1], dtype=float)
        delta_r2 = full_baseline_r2 - leave_out_r2

        for site_idx, record in enumerate(records):
            rows.append(
                {
                    "site_id": record.site_id,
                    "name": record.name,
                    "site_class": record.site_class,
                    "basin": record.basin,
                    "basis_name": basis_name,
                    "basis_coefficient": float(basis_coeffs[basis_idx, site_idx]),
                    "basis_abs_coefficient": float(abs(basis_coeffs[basis_idx, site_idx])),
                    "basis_single_term_model_r2_vs_raw": float(single_r2_raw[site_idx]),
                    "basis_single_term_model_coefficient": (
                        float(single_coeffs[-1, site_idx]) if basis_idx != 0 else float(single_coeffs[0, site_idx])
                    ),
                    "basis_drop_one_delta_r2_vs_raw": float(delta_r2[site_idx]),
                }
            )

        summary_rows.append(
            {
                "basis_name": basis_name,
                "mean_abs_coefficient": float(np.mean(np.abs(basis_coeffs[basis_idx]))),
                "median_abs_coefficient": float(np.median(np.abs(basis_coeffs[basis_idx]))),
                "mean_single_term_model_r2_vs_raw": float(np.mean(single_r2_raw)),
                "median_single_term_model_r2_vs_raw": float(np.median(single_r2_raw)),
                "mean_drop_one_delta_r2_vs_raw": float(np.mean(delta_r2)),
                "median_drop_one_delta_r2_vs_raw": float(np.median(delta_r2)),
            }
        )

    detail_df = pd.DataFrame(rows).sort_values(["basis_name", "site_class", "site_id"])
    summary_df = pd.DataFrame(summary_rows).sort_values("mean_drop_one_delta_r2_vs_raw", ascending=False)
    return detail_df, summary_df


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    config_path = shared.resolve_user_path(args.config)
    output_dir = shared.resolve_user_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_config(config_path, args)
    base_dir = (repo_root / cfg["data"]["base_dir"]).resolve()
    id_file = (repo_root / cfg["data"]["id_file"]).resolve()

    id_map = shared.load_id_map(id_file)
    records = shared.load_site_records(base_dir, id_map, cfg)
    shared.write_metadata(records, output_dir)

    times, raw_matrix = shared.build_aligned_matrix(records)
    mean_cfg = cfg["data"].get("mean_forcing_hint", {})
    mean_signal_raw, mean_signal_norm, hint_info = shared.align_hint_series(times, mean_cfg, repo_root)
    mean_harmonic_factors = shared.get_mean_forcing_harmonic_factors(mean_cfg)
    mean_basis_names = shared.get_mean_forcing_basis_names(mean_harmonic_factors)
    site_weights = shared.build_site_weights(records, float(cfg["training"].get("record_length_weight_power", 0.0)))
    annual_latent_dim = int(cfg["annual_model"]["latent_dim"])

    mean_baseline_full = np.zeros_like(raw_matrix)
    baseline_coeffs = None
    if hint_info.get("enabled"):
        ridge = float(mean_cfg.get("composite_ridge", 1.0e-4))
        l1_alpha = float(mean_cfg.get("lasso_alpha", 0.0))
        fit_indices = np.arange(len(times), dtype=int)
        mean_baseline_full, baseline_coeffs = shared.fit_mean_forcing_baseline(
            mean_signal_raw,
            mean_signal_norm,
            raw_matrix,
            fit_indices,
            ridge=ridge,
            l1_alpha=l1_alpha,
            harmonic_factors=mean_harmonic_factors,
        )

    residual_matrix = raw_matrix - mean_baseline_full
    annual_years, annual_levels, months_per_year, year_inverse = build_annual_level_matrix(times, residual_matrix)
    annual_fit = fit_annual_increment_factor_model(annual_levels, site_weights, annual_latent_dim)

    annual_component_full = expand_annual_levels(annual_fit["annual_level_recon"], year_inverse)
    recon_full = mean_baseline_full + annual_component_full
    residual_truth = residual_matrix
    monthly_latent_levels = annual_fit["latent_levels"][year_inverse]

    metrics_rows = []
    final_rmse, final_mae, final_corr = shared.summarize_metrics(raw_matrix, recon_full)
    baseline_rmse, baseline_mae, baseline_corr = shared.summarize_metrics(raw_matrix, mean_baseline_full)
    annual_rmse, annual_mae, annual_corr = shared.summarize_metrics(residual_truth, annual_component_full)
    final_r2 = coefficient_of_determination(raw_matrix, recon_full)
    baseline_r2 = coefficient_of_determination(raw_matrix, mean_baseline_full)
    annual_r2_on_residual = coefficient_of_determination(residual_truth, annual_component_full)
    annual_r2_on_raw = coefficient_of_determination(raw_matrix, annual_component_full)
    annual_corr_increment = final_corr - baseline_corr
    baseline_corr_fraction = shared.fractional_corr_contribution(baseline_corr, final_corr)
    annual_corr_fraction = shared.fractional_corr_contribution(annual_corr_increment, final_corr)
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
                "r2": float(final_r2[idx]),
                "baseline_rmse": float(baseline_rmse[idx]),
                "baseline_mae": float(baseline_mae[idx]),
                "baseline_corr": float(baseline_corr[idx]),
                "baseline_r2": float(baseline_r2[idx]),
                "baseline_corr_fraction_of_final": float(baseline_corr_fraction[idx]),
                "annual_component_rmse": float(annual_rmse[idx]),
                "annual_component_mae": float(annual_mae[idx]),
                "annual_component_corr": float(annual_corr[idx]),
                "annual_component_r2_on_residual": float(annual_r2_on_residual[idx]),
                "annual_component_r2_on_raw": float(annual_r2_on_raw[idx]),
                "annual_corr_increment": float(annual_corr_increment[idx]),
                "annual_corr_fraction_of_final": float(annual_corr_fraction[idx]),
                "record_months": int(record.record_months),
                "record_years": float(record.record_years),
                "training_weight": float(site_weights[idx]),
            }
        )
    metrics_df = pd.DataFrame(metrics_rows).sort_values(["site_class", "corr"], ascending=[True, False])
    metrics_df.to_csv(output_dir / "reconstruction_metrics.csv", index=False)

    basis_detail_df = None
    basis_summary_df = None
    if baseline_coeffs is not None and mean_signal_norm is not None:
        basis_all = shared.build_mean_forcing_basis(mean_signal_norm, harmonic_factors=mean_harmonic_factors)
        basis_detail_df, basis_summary_df = build_mean_forcing_basis_determination(
            basis_all,
            baseline_coeffs,
            raw_matrix,
            mean_baseline_full,
            records,
            mean_basis_names,
            float(mean_cfg.get("composite_ridge", 1.0e-4)),
            float(mean_cfg.get("lasso_alpha", 0.0)),
        )
        basis_detail_df.to_csv(output_dir / "mean_forcing_basis_determination.csv", index=False)
        basis_summary_df.to_csv(output_dir / "mean_forcing_basis_determination_summary.csv", index=False)

    latent_df = pd.DataFrame({"time": times})
    for dim in range(monthly_latent_levels.shape[1]):
        latent_df[f"z{dim + 1}"] = monthly_latent_levels[:, dim]
    latent_df.to_csv(output_dir / "latent_states.csv", index=False)
    write_annual_latent_states(
        output_dir / "annual_latent_states.csv",
        annual_years,
        annual_fit["latent_levels"],
        annual_fit["latent_increments"],
    )

    latent_mean_df, latent_mean_aligned_df = shared.compare_latent_to_mean_forcing(
        times,
        mean_signal_raw,
        mean_signal_norm,
        np.arange(len(times), dtype=int),
        monthly_latent_levels,
    )
    if latent_mean_df is not None and latent_mean_aligned_df is not None:
        latent_mean_df.to_csv(output_dir / "latent_mean_forcing_comparison.csv", index=False)
        latent_mean_aligned_df.to_csv(output_dir / "latent_mean_forcing_aligned.csv", index=False)

    write_wide_matrix(output_dir / "aligned_model_matrix.csv", times, records, raw_matrix)
    write_wide_matrix(output_dir / "annual_component_matrix.csv", times, records, annual_component_full)
    if hint_info.get("enabled"):
        pd.DataFrame({"time": times, "mean_forcing_raw": mean_signal_raw, "mean_forcing_normalized": mean_signal_norm}).to_csv(
            output_dir / "aligned_mean_forcing_hint.csv",
            index=False,
        )
        write_wide_matrix(output_dir / "mean_forcing_baseline_matrix.csv", times, records, mean_baseline_full)

    annual_levels_df = pd.DataFrame({"year": annual_years, "months_in_year": months_per_year})
    for idx, record in enumerate(records):
        annual_levels_df[record.site_id] = annual_fit["annual_level_recon"][:, idx]
    annual_levels_df.to_csv(output_dir / "annual_level_matrix.csv", index=False)

    annual_increment_df = pd.DataFrame({"year": annual_years[1:]})
    for idx, record in enumerate(records):
        annual_increment_df[record.site_id] = annual_fit["annual_increment_recon"][:, idx]
    annual_increment_df.to_csv(output_dir / "annual_increment_matrix.csv", index=False)

    write_site_loadings(output_dir / "annual_site_loadings.csv", records, annual_fit["site_loadings"])

    history_df = pd.DataFrame(
        {
            "component": np.arange(1, len(annual_fit["singular_values"]) + 1, dtype=int),
            "singular_value": annual_fit["singular_values"],
            "explained_variance_ratio": annual_fit["explained_variance_ratio"],
            "cumulative_explained_variance_ratio": np.cumsum(annual_fit["explained_variance_ratio"]),
        }
    )
    history_df.to_csv(output_dir / "training_history.csv", index=False)

    peaks = shared.top_spectral_peaks(
        monthly_latent_levels,
        dt_years=1.0 / 12.0,
        top_n=int(cfg["report"]["top_spectral_peaks"]),
    )
    peaks = shared.attach_period_matches(peaks, cfg["report"]["period_library_years"])
    peaks.to_csv(output_dir / "latent_spectral_peaks.csv", index=False)
    shared.spectral_entropy(monthly_latent_levels).to_csv(output_dir / "latent_spectral_entropy.csv", index=False)

    shared.plot_latent_spectra(monthly_latent_levels, dt_years=1.0 / 12.0, output_path=output_dir / "latent_spectra.png")
    shared.plot_reconstruction_examples(
        times,
        raw_matrix,
        recon_full,
        records,
        metrics_df,
        output_dir / "reconstruction_examples.png",
        per_class=int(cfg["report"]["example_sites_per_class"]),
    )
    shared.plot_latent_mean_forcing_comparison(
        latent_mean_aligned_df,
        latent_mean_df,
        output_dir / "latent_mean_forcing_comparison.png",
    )

    mean_final_corr = float(metrics_df["corr"].mean())
    mean_baseline_corr = float(metrics_df["baseline_corr"].mean())
    mean_annual_corr = float(metrics_df["annual_component_corr"].mean())
    mean_annual_increment = mean_final_corr - mean_baseline_corr
    mean_final_r2 = float(metrics_df["r2"].mean())
    mean_baseline_r2 = float(metrics_df["baseline_r2"].mean())
    mean_annual_r2_on_residual = float(metrics_df["annual_component_r2_on_residual"].mean())
    mean_annual_r2_on_raw = float(metrics_df["annual_component_r2_on_raw"].mean())
    summary = {
        "config_path": str(config_path),
        "output_dir": str(output_dir),
        "approach": "annual_increment_factor",
        "site_count": len(records),
        "site_classes": dict(Counter(record.site_class for record in records)),
        "excluded_prefixes": cfg["data"].get("exclude_prefixes", []),
        "excluded_sites": cfg["data"].get("exclude_sites", []),
        "mean_forcing_hint": hint_info,
        "record_length_weight_power": float(cfg["training"].get("record_length_weight_power", 0.0)),
        "shared_span": {"start_year": float(times[0]), "stop_year": float(times[-1]), "months": int(len(times))},
        "annual_year_count": int(len(annual_years)),
        "annual_months_per_year_min": int(months_per_year.min()),
        "annual_months_per_year_max": int(months_per_year.max()),
        "annual_latent_dim": int(monthly_latent_levels.shape[1]),
        "mean_reconstruction_corr": mean_final_corr,
        "median_reconstruction_corr": float(metrics_df["corr"].median()),
        "mean_reconstruction_r2": mean_final_r2,
        "median_reconstruction_r2": float(metrics_df["r2"].median()),
        "mean_baseline_corr": mean_baseline_corr,
        "median_baseline_corr": float(metrics_df["baseline_corr"].median()),
        "mean_baseline_r2": mean_baseline_r2,
        "median_baseline_r2": float(metrics_df["baseline_r2"].median()),
        "mean_annual_component_corr": mean_annual_corr,
        "median_annual_component_corr": float(metrics_df["annual_component_corr"].median()),
        "mean_annual_component_r2_on_residual": mean_annual_r2_on_residual,
        "median_annual_component_r2_on_residual": float(metrics_df["annual_component_r2_on_residual"].median()),
        "mean_annual_component_r2_on_raw": mean_annual_r2_on_raw,
        "median_annual_component_r2_on_raw": float(metrics_df["annual_component_r2_on_raw"].median()),
        "period_library_years": cfg["report"]["period_library_years"],
        "corr_contribution_breakdown": {
            "combined_mean_corr": mean_final_corr,
            "mean_forcing_mean_corr": mean_baseline_corr,
            "annual_component_mean_corr": mean_annual_corr,
            "annual_increment_mean_corr": mean_annual_increment,
            "mean_forcing_fraction_of_combined_corr": (
                float(mean_baseline_corr / mean_final_corr) if abs(mean_final_corr) > 1.0e-8 else None
            ),
            "annual_increment_fraction_of_combined_corr": (
                float(mean_annual_increment / mean_final_corr) if abs(mean_final_corr) > 1.0e-8 else None
            ),
        },
        "r2_component_breakdown": {
            "combined_mean_r2": mean_final_r2,
            "mean_forcing_mean_r2": mean_baseline_r2,
            "annual_component_mean_r2_on_residual": mean_annual_r2_on_residual,
            "annual_component_mean_r2_on_raw": mean_annual_r2_on_raw,
        },
    }
    if baseline_coeffs is not None:
        summary["mean_forcing_basis_coeff_shape"] = list(np.asarray(baseline_coeffs).shape)
    if basis_summary_df is not None and not basis_summary_df.empty:
        summary["mean_forcing_basis_determination_summary"] = {
            row["basis_name"]: {
                "mean_single_term_model_r2_vs_raw": float(row["mean_single_term_model_r2_vs_raw"]),
                "mean_drop_one_delta_r2_vs_raw": float(row["mean_drop_one_delta_r2_vs_raw"]),
            }
            for _, row in basis_summary_df.iterrows()
        }
    if not history_df.empty:
        summary["annual_cumulative_explained_variance_ratio"] = float(history_df["cumulative_explained_variance_ratio"].iloc[-1])
    if latent_mean_df is not None and not latent_mean_df.empty:
        best_latent_row = latent_mean_df.iloc[0]
        summary["best_latent_mean_forcing_match"] = {
            "latent_dim": int(best_latent_row["latent_dim"]),
            "corr_with_mean_forcing": float(best_latent_row["corr_with_mean_forcing"]),
            "abs_corr_with_mean_forcing": float(best_latent_row["abs_corr_with_mean_forcing"]),
        }
    (output_dir / "run_summary.json").write_text(json.dumps(shared.json_ready(summary), indent=2), encoding="utf-8")

    print(f"Wrote artifacts to {output_dir}")
    print(f"Sites used: {len(records)}")
    print(f"Shared span: {times[0]:.3f} -> {times[-1]:.3f} ({len(times)} months)")
    print(f"Annual years modeled: {annual_years[0]} -> {annual_years[-1]} ({len(annual_years)} years)")
    print(f"Mean reconstruction correlation: {mean_final_corr:.4f}")
    print(
        "Mean CC breakdown: "
        f"mean_forcing.dat={mean_baseline_corr:.4f} "
        f"({summary['corr_contribution_breakdown']['mean_forcing_fraction_of_combined_corr']:.1%} of combined), "
        f"annual_increment={mean_annual_increment:.4f} "
        f"({summary['corr_contribution_breakdown']['annual_increment_fraction_of_combined_corr']:.1%} of combined)"
    )


if __name__ == "__main__":
    main()
