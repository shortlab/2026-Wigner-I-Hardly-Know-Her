from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io_parser import parse_segments
from ..methods import find_lefts_rights_mc, hook_mc_qad, spline_integral_orig


def _heating_rate(seg) -> float:
    t = seg.time_s
    tr = seg.tr_c
    if len(t) > 25 and np.isfinite(t[24]) and np.isfinite(t[23]) and (t[24] != t[23]):
        return float((tr[24] - tr[23]) / (t[24] - t[23]))
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
) -> dict:
    fp = Path(file_path)
    segs = parse_segments(fp)
    if not segs:
        raise RuntimeError(f"No segments parsed from {fp}")

    heat = [s for s in segs if s.direction == "heating"]
    use = heat if heat else segs
    if len(use) == 0:
        raise RuntimeError("No usable segments")

    out_dir = Path(output_dir)
    if not out_dir.is_absolute():
        out_dir = fp.parent / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    h_fusion = 28.6
    heat_capacity = 0.23

    for i, seg in enumerate(use):
        hr = _heating_rate(seg)
        bounds = find_lefts_rights_mc(seg.tr_c, seg.q_mw, n_samples=n_mc, seed=seed + i)
        if bounds is None:
            continue

        hcrys, _, _ = spline_integral_orig(seg.tr_c, seg.q_mw, bounds["ll"], bounds["lr"], bounds["rl"], bounds["rr"])
        hcrys_mean = float(np.nanmean(hcrys))
        hcrys_std = float(np.nanstd(hcrys))

        t_lr_mean = float(np.nanmean(bounds["lr"]))
        excl = seg.tr_c < t_lr_mean
        hookh = np.nan
        if np.count_nonzero(excl) >= 10:
            h_start, h_end, _, _, _ = hook_mc_qad(seg.tr_c[excl], seg.q_mw[excl], fig_on=False, n=n_mc, seed=seed + i)
            if np.isfinite(h_start[0]) and np.isfinite(h_end[0]):
                hookh = float(abs(h_end[0] - h_start[0]))

        mass_hcrys = float(abs(hcrys_mean * 1e-3 / hr / h_fusion / 1e-9)) if np.isfinite(hr) and hr != 0 else np.nan
        mass_hook = float(abs(hookh * 1e-3 / hr / heat_capacity / 1e-9)) if np.isfinite(hr) and hr != 0 and np.isfinite(hookh) else np.nan

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
            }
        )

        plt.figure(figsize=(8, 4))
        plt.plot(seg.tr_c, seg.q_mw, linewidth=0.8)
        plt.axvline(bounds["lrcenter"], linestyle="--", linewidth=0.8)
        plt.axvline(bounds["rlcenter"], linestyle="--", linewidth=0.8)
        plt.title(f"Segment {i} - seg_id {seg.seg_id}")
        plt.xlabel("Temperature (C)")
        plt.ylabel("Power (mW)")
        plt.tight_layout()
        plt.savefig(out_dir / f"segment_{i:02d}.png", dpi=150)
        plt.close()

    if not records:
        raise RuntimeError("No successful analysis records generated")

    df = pd.DataFrame(records)
    csv_path = out_dir / "indium_analysis_summary.csv"
    df.to_csv(csv_path, index=False)

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
    }
