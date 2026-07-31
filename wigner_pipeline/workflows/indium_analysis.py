from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io_parser import parse_segments
from ..methods import find_lefts_rights_mc, hook_mc_qad, spline_integral_orig


def _heating_rate(seg) -> float:
    """Estimate segment heating rate in K/s with MATLAB-parity fallback behavior."""
    # Prefer the same local finite-difference style used in MATLAB (around index ~25)
    # to keep parity with legacy heating-rate estimates.
    t = seg.time_s
    tr = seg.tr_c
    if len(t) > 25 and np.isfinite(t[24]) and np.isfinite(t[23]) and (t[24] != t[23]):
        return float((tr[24] - tr[23]) / (t[24] - t[23]))
    # Fallback for short/irregular segments: robust median gradient.
    g = np.gradient(tr, t)
    g = g[np.isfinite(g)]
    if len(g) == 0:
        return float("nan")
    return float(np.nanmedian(g))


def run_indium_analysis(
    file_path: str | Path,
    output_dir: str | Path = "outputs/indium_analysis",
    n_mc: int = 1000,
    seed: int = 42,
    fixed_lr_center: float | None = None,
    fixed_rl_center: float | None = None,
    fixed_outer_pad_c: float = 15.0,
    baseline_fit_offset_mw: float = 1.0,
) -> dict:
    """Run the indium mass-estimation workflow and persist summary/diagnostic outputs."""
    fp = Path(file_path)

    use_fixed_bounds = (fixed_lr_center is not None) or (fixed_rl_center is not None)
    if use_fixed_bounds and (fixed_lr_center is None or fixed_rl_center is None):
        raise ValueError("Both fixed_lr_center and fixed_rl_center must be provided together")
    if use_fixed_bounds and fixed_lr_center >= fixed_rl_center:
        raise ValueError("fixed_lr_center must be less than fixed_rl_center")

    # Parse all segments from one exported text file.
    segs = parse_segments(fp)
    if not segs:
        raise RuntimeError(f"No segments parsed from {fp}")

    # Prefer heating segments for indium melting analysis.
    heat = [s for s in segs if s.direction == "heating"]
    use = heat if heat else segs
    if len(use) == 0:
        raise RuntimeError("No usable segments")

    out_dir = Path(output_dir)
    if not out_dir.is_absolute():
        out_dir = fp.parent / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    # Constants carried from MATLAB workflow.
    # Reference values for indium used in legacy script:
    #   Enthalpy of fusion: 3283 J/mol (28.6 J/g)
    #   Cp at 300K: 26.98 J/K/mol (0.23 J/K/g)
    # Source: https://akjournals.com/view/journals/10973/13/3/article-p419.xml
    h_fusion = 28.6
    heat_capacity = 0.23

    for i, seg in enumerate(use):
        # Per-segment heating rate (K/s).
        hr = _heating_rate(seg)

        # Detect MC integration bounds around the transition, or use fixed
        # user-provided centers for sensitivity tests.
        if use_fixed_bounds:
            lr_center = float(fixed_lr_center)
            rl_center = float(fixed_rl_center)
            bounds = {
                "ll": np.full(n_mc, lr_center - float(fixed_outer_pad_c), dtype=float),
                "lr": np.full(n_mc, lr_center, dtype=float),
                "rl": np.full(n_mc, rl_center, dtype=float),
                "rr": np.full(n_mc, rl_center + float(fixed_outer_pad_c), dtype=float),
                "dTmax": np.nan,
                "dTmin": np.nan,
                "lrcenter": lr_center,
                "rlcenter": rl_center,
            }
            bounds_mode = "fixed"
        else:
            bounds = find_lefts_rights_mc(seg.tr_c, seg.q_mw, n_samples=n_mc, seed=seed + i)
            if bounds is None:
                continue
            bounds_mode = "auto"

        # Enthalpy-like area from baseline-corrected integration across MC draws.
        hcrys, _, _ = spline_integral_orig(
            seg.tr_c,
            seg.q_mw,
            bounds["ll"],
            bounds["lr"],
            bounds["rl"],
            bounds["rr"],
            baseline_fit_offset_mw=baseline_fit_offset_mw,
        )
        hcrys_mean = float(np.nanmean(hcrys))
        hcrys_std = float(np.nanstd(hcrys))

        # Hook method is computed on pre-transition data (left of lr center).
        t_lr_mean = float(np.nanmean(bounds["lr"]))
        excl = seg.tr_c < t_lr_mean
        hookh = np.nan
        hook_start_xy = (np.nan, np.nan)
        hook_end_xy = (np.nan, np.nan)
        if np.count_nonzero(excl) >= 10:
            t_hook = np.asarray(seg.tr_c[excl], dtype=float)
            q_hook = np.asarray(seg.q_mw[excl], dtype=float)
            h_start, h_end, detail, _, _ = hook_mc_qad(t_hook, q_hook, fig_on=False, n=n_mc, seed=seed + i)
            if np.isfinite(h_start[0]) and np.isfinite(h_end[0]):
                hookh = float(abs(h_end[0] - h_start[0]))

            # Diagnostic marker locations for visualization on the segment plot.
            mode = str(detail[0]) if detail and detail[0] is not None else "heat"
            x_intersect = float(detail[2][0]) if detail and len(detail) >= 3 and np.ndim(detail[2]) > 0 else np.nan
            if np.isfinite(h_start[0]) and len(t_hook) > 0:
                x_start = float(t_hook[0]) if mode == "heat" else float(t_hook[-1])
                hook_start_xy = (x_start, float(h_start[0]))
            if np.isfinite(h_end[0]) and np.isfinite(x_intersect):
                hook_end_xy = (x_intersect, float(h_end[0]))

        # Convert each method-specific signal into mass estimate (ng).
        # hcrys_mean is mW*C integrated and then scaled by heating rate and latent heat.
        mass_hcrys = float(abs(hcrys_mean * 1e-3 / hr / h_fusion / 1e-9)) if np.isfinite(hr) and hr != 0 else np.nan
        # hookh is mW offset scaled by heating rate and heat capacity.
        mass_hook = float(abs(hookh * 1e-3 / hr / heat_capacity / 1e-9)) if np.isfinite(hr) and hr != 0 and np.isfinite(hookh) else np.nan

        # Keep compact per-segment summary for CSV and downstream comparisons.
        records.append(
            {
                "segment": i,
                "seg_id": seg.seg_id,
                "direction": seg.direction,
                "n_points": len(seg.tr_c),
                "heating_rate_K_per_s": hr,
                "hcrys_mean": hcrys_mean,
                "hcrys_std": hcrys_std,
                "hook_height": hookh,
                "mass_hcrys_ng": mass_hcrys,
                "mass_hook_ng": mass_hook,
                "lr_center_C": bounds["lrcenter"],
                "rl_center_C": bounds["rlcenter"],
                "bounds_mode": bounds_mode,
                "baseline_fit_offset_mW": float(baseline_fit_offset_mw),
            }
        )

        # Per-segment diagnostic plot showing raw trace, detected bounds, and hook markers.
        plt.figure(figsize=(8, 4))
        plt.plot(seg.tr_c, seg.q_mw, linewidth=0.8)
        plt.axvline(bounds["lrcenter"], linestyle="--", linewidth=0.8)
        plt.axvline(bounds["rlcenter"], linestyle="--", linewidth=0.8)

        if np.isfinite(hook_start_xy[0]) and np.isfinite(hook_start_xy[1]):
            plt.scatter(
                [hook_start_xy[0]],
                [hook_start_xy[1]],
                s=20,
                marker="o",
                label="hook start",
                zorder=4,
            )
            plt.annotate(
                f"start ({hook_start_xy[0]:.2f}, {hook_start_xy[1]:.4f})",
                xy=(hook_start_xy[0], hook_start_xy[1]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=7,
            )
        if np.isfinite(hook_end_xy[0]) and np.isfinite(hook_end_xy[1]):
            plt.scatter(
                [hook_end_xy[0]],
                [hook_end_xy[1]],
                s=24,
                marker="x",
                linewidths=1.0,
                label="hook end",
                zorder=5,
            )
            plt.annotate(
                f"end ({hook_end_xy[0]:.2f}, {hook_end_xy[1]:.4f})",
                xy=(hook_end_xy[0], hook_end_xy[1]),
                xytext=(6, -10),
                textcoords="offset points",
                fontsize=7,
            )

        plt.title(f"Segment {i} - seg_id {seg.seg_id}")
        plt.xlabel("Temperature (C)")
        plt.ylabel("Power (mW)")
        plt.legend(loc="best", fontsize=7)
        plt.tight_layout()
        plt.savefig(out_dir / f"segment_{i:02d}.png", dpi=150)
        plt.close()

    if not records:
        raise RuntimeError("No successful analysis records generated")

    # Persist numerical outputs for auditing and cross-run comparisons.
    df = pd.DataFrame(records)
    csv_path = out_dir / "indium_analysis_summary.csv"
    df.to_csv(csv_path, index=False)

    # Plot only records with finite estimates and valid positive heating rates.
    valid = df[np.isfinite(df["mass_hcrys_ng"]) & np.isfinite(df["mass_hook_ng"]) & np.isfinite(df["heating_rate_K_per_s"]) & (df["heating_rate_K_per_s"] > 0)]
    if not valid.empty:
        plt.figure(figsize=(8, 5))
        plt.semilogx(valid["heating_rate_K_per_s"], valid["mass_hcrys_ng"], ".", markersize=10, label="enthalpy")
        plt.semilogx(valid["heating_rate_K_per_s"], valid["mass_hook_ng"], "x", markersize=8, label="hook")
        plt.xlabel("Heating rate (K/s)")
        plt.ylabel("Estimated mass (ng)")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(out_dir / "mass_vs_heating_rate.png", dpi=150)
        plt.close()

    return {
        "file": str(fp),
        "n_segments_input": len(segs),
        "n_segments_used": len(use),
        "n_records": len(df),
        "output_dir": str(out_dir),
        "csv": str(csv_path),
        "bounds_mode": "fixed" if use_fixed_bounds else "auto",
        "fixed_lr_center_C": float(fixed_lr_center) if use_fixed_bounds else None,
        "fixed_rl_center_C": float(fixed_rl_center) if use_fixed_bounds else None,
        "fixed_outer_pad_C": float(fixed_outer_pad_c) if use_fixed_bounds else None,
        "baseline_fit_offset_mW": float(baseline_fit_offset_mw),
    }
