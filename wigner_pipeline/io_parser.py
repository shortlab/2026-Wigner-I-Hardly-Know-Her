from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .models import Segment


_CURVE_VALUES_RE = re.compile(r"Curve Values:", re.IGNORECASE)
_SEG_ID_RE = re.compile(r"\$\]\s*(\d+)\s*\[")
_MASS_RE = re.compile(r",\s*([\d.]+)\s*mg", re.IGNORECASE)


def _to_float(token: str) -> float:
    """Parse numeric tokens while preserving explicit NaN markers."""
    token = token.strip()
    if token.lower() == "nan":
        return float("nan")
    return float(token)


def _is_data_line(line: str) -> bool:
    """Return True when a text line starts with an integer row index."""
    s = line.strip()
    if not s:
        return False
    first = s.split()[0]
    try:
        int(first)
        return True
    except ValueError:
        return False


def _parse_stare_blocks(path: Path, lines: list[str]) -> list[Segment]:
    """Parse STAR(e)-style exported blocks delimited by 'Curve Values:' headers."""
    starts = [i for i, ln in enumerate(lines) if _CURVE_VALUES_RE.search(ln)]
    if not starts:
        return []

    segments: list[Segment] = []
    for file_order, start in enumerate(starts):
        seg_id = file_order
        mass_mg = None
        for back in range(1, 12):
            idx = start - back
            if idx < 0:
                break
            mseg = _SEG_ID_RE.search(lines[idx])
            if mseg:
                seg_id = int(mseg.group(1))
            mmass = _MASS_RE.search(lines[idx])
            if mmass and mass_mg is None:
                mass_mg = float(mmass.group(1))

        rows: list[list[float]] = []
        for ln in lines[start + 1 :]:
            if not _is_data_line(ln):
                if rows:
                    break
                continue
            parts = ln.split()
            vals = []
            for p in parts:
                try:
                    vals.append(_to_float(p))
                except ValueError:
                    break
            if len(vals) >= 4:
                rows.append(vals)

        if not rows:
            continue

        arr = np.asarray(rows, dtype=float)
        n_cols = arr.shape[1]
        if n_cols >= 5:
            time_s = arr[:, 1]
            ts_c = arr[:, n_cols - 3]
            tr_c = arr[:, n_cols - 2]
            q_mw = arr[:, n_cols - 1]
        else:
            time_s = np.arange(arr.shape[0], dtype=float)
            ts_c = arr[:, n_cols - 3]
            tr_c = arr[:, n_cols - 2]
            q_mw = arr[:, n_cols - 1]

        valid = np.isfinite(ts_c) & np.isfinite(tr_c) & np.isfinite(q_mw)
        time_s = np.asarray(time_s[valid], dtype=float)
        ts_c = np.asarray(ts_c[valid], dtype=float)
        tr_c = np.asarray(tr_c[valid], dtype=float)
        q_mw = np.asarray(q_mw[valid], dtype=float)

        if len(time_s) < 2:
            continue

        direction = "cooling" if tr_c[0] > tr_c[-1] else "heating"
        segments.append(
            Segment(
                source_file=path,
                segment_index=file_order,
                seg_id=seg_id,
                direction=direction,
                time_s=time_s,
                ts_c=ts_c,
                tr_c=tr_c,
                q_mw=q_mw,
                mass_mg=mass_mg,
            )
        )

    segments.sort(key=lambda x: x.seg_id)
    return segments


def _parse_multipart_numeric(path: Path, lines: list[str]) -> list[Segment]:
    """Parse plain numeric exports that separate segments by index reset to 0."""
    rows: list[list[float]] = []
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        parts = s.split()
        vals = []
        for p in parts:
            try:
                vals.append(_to_float(p))
            except ValueError:
                vals = []
                break
        if len(vals) >= 4:
            rows.append(vals)

    if not rows:
        return []

    segments: list[list[list[float]]] = []
    current: list[list[float]] = []
    for row in rows:
        idx = row[0]
        if np.isnan(idx):
            continue
        if idx == 0 and current:
            segments.append(current)
            current = []
        current.append(row)
    if current:
        segments.append(current)

    out: list[Segment] = []
    for i, seg_rows in enumerate(segments):
        arr = np.asarray(seg_rows, dtype=float)
        n_cols = arr.shape[1]
        if n_cols >= 5:
            time_s = arr[:, 1]
            ts_c = arr[:, n_cols - 3]
            tr_c = arr[:, n_cols - 2]
            q_mw = arr[:, n_cols - 1]
        else:
            time_s = arr[:, 1] if n_cols >= 4 else np.arange(arr.shape[0], dtype=float)
            ts_c = arr[:, n_cols - 3]
            tr_c = arr[:, n_cols - 2]
            q_mw = arr[:, n_cols - 1]

        valid = np.isfinite(ts_c) & np.isfinite(tr_c) & np.isfinite(q_mw)
        time_s = np.asarray(time_s[valid], dtype=float)
        ts_c = np.asarray(ts_c[valid], dtype=float)
        tr_c = np.asarray(tr_c[valid], dtype=float)
        q_mw = np.asarray(q_mw[valid], dtype=float)
        if len(time_s) < 2:
            continue
        direction = "cooling" if tr_c[0] > tr_c[-1] else "heating"
        out.append(
            Segment(
                source_file=path,
                segment_index=i,
                seg_id=i,
                direction=direction,
                time_s=time_s,
                ts_c=ts_c,
                tr_c=tr_c,
                q_mw=q_mw,
                mass_mg=None,
            )
        )
    return out


def parse_segments(path: str | Path) -> list[Segment]:
    """Parse one text export into ordered heating/cooling Segment objects."""
    p = Path(path)
    lines = p.read_text(encoding="latin-1", errors="replace").splitlines()
    segments = _parse_stare_blocks(p, lines)
    if segments:
        return segments
    return _parse_multipart_numeric(p, lines)


def parser_smoke(path: str | Path) -> dict:
    """Return a compact parse summary for quick parser validation/debugging."""
    segs = parse_segments(path)
    if not segs:
        raise RuntimeError(f"No segments parsed from {path}")
    return {
        "file": str(path),
        "n_segments": len(segs),
        "segment_ids": [s.seg_id for s in segs],
        "directions": [s.direction for s in segs],
        "lengths": [int(len(s.tr_c)) for s in segs],
    }
