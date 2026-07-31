from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from ..io_parser import parse_segments


def _stack_file_segments(file_segments: list) -> tuple[np.ndarray, np.ndarray]:
    if len(file_segments) < 1:
        raise RuntimeError("No segments in file")
    min_len = min(len(s.ts_c) for s in file_segments)
    t = np.column_stack([s.ts_c[:min_len] for s in file_segments])
    q = np.column_stack([s.q_mw[:min_len] for s in file_segments])
    return t, q


def run_interface_experiment(
    workspace_dir: str | Path,
    baseline_pattern: str = "cle_ind_baseline",
    output_dir: str | Path = "outputs/interface_experiment",
) -> dict:
    root = Path(workspace_dir)
    files = sorted(root.glob("*.txt"))
    if not files:
        raise RuntimeError("No .txt files found")

    output_path = Path(output_dir)
    if not output_path.is_absolute():
        output_path = root / output_path
    output_path.mkdir(parents=True, exist_ok=True)

    records = []
    for fp in files:
        segs = parse_segments(fp)
        if not segs:
            continue
        temp_mat, flow_mat = _stack_file_segments(segs)
        records.append({"file": fp, "temp": temp_mat, "flow": flow_mat, "n_segments": flow_mat.shape[1]})

    if not records:
        raise RuntimeError("No parseable files found")

    baseline_matches = [r for r in records if baseline_pattern.lower() in r["file"].name.lower()]
    if len(baseline_matches) != 1:
        raise RuntimeError(f"Expected exactly one baseline file containing '{baseline_pattern}', found {len(baseline_matches)}")
    baseline = baseline_matches[0]

    global_min, global_max = np.inf, -np.inf
    plotted = []
    for rec in records:
        is_ind_file = ("ind" in rec["file"].name.lower()) and (rec["file"] != baseline["file"])
        seg_y = []
        for s in range(rec["n_segments"]):
            x = rec["temp"][:, s]
            y = rec["flow"][:, s]
            if is_ind_file and s < baseline["n_segments"]:
                f = interp1d(baseline["temp"][:, s], baseline["flow"][:, s], kind="linear", bounds_error=False, fill_value=np.nan)
                y = y - f(x)
            seg_y.append(y)
            global_min = min(global_min, float(np.nanmin(y)))
            global_max = max(global_max, float(np.nanmax(y)))
        plotted.append({"rec": rec, "is_ind": is_ind_file, "seg_y": seg_y})

    ypad = 0.05 * (global_max - global_min if global_max > global_min else 1.0)
    ylims = (global_min - ypad, global_max + ypad)

    output_files = []
    for i, item in enumerate(plotted, start=1):
        rec = item["rec"]
        plt.figure(figsize=(8, 5))
        for s in range(rec["n_segments"]):
            plt.plot(rec["temp"][:, s], item["seg_y"][s], linewidth=0.8, label=f"seg_{s+1}")
        plt.xlim(20, 325)
        plt.ylim(*ylims)
        plt.xlabel("Temperature (C)")
        plt.ylabel("Delta Power (mW)" if item["is_ind"] else "Power (mW)")
        title = rec["file"].name + (" (indium baseline subtracted)" if item["is_ind"] else "")
        plt.title(title)
        if item["is_ind"]:
            plt.axhline(0.0, linestyle="--", linewidth=0.8)
        plt.legend(loc="best", fontsize=7)
        plt.tight_layout()
        out = output_path / f"{i:02d}_{rec['file'].stem}.png"
        plt.savefig(out, dpi=150)
        plt.close()
        output_files.append(str(out))

    return {
        "n_files": len(records),
        "baseline_file": baseline["file"].name,
        "output_dir": str(output_path),
        "plots": output_files,
    }
