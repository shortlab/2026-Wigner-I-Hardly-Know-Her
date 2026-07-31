from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Segment:
    """One parsed DSC segment with metadata and aligned numeric vectors."""
    source_file: Path
    segment_index: int
    seg_id: int
    direction: str
    time_s: np.ndarray
    ts_c: np.ndarray
    tr_c: np.ndarray
    q_mw: np.ndarray
    mass_mg: float | None = None

    @property
    def label(self) -> str:
        """Human-readable identifier used in logs/plots."""
        return f"{self.source_file.stem} [seg_id {self.seg_id} - {self.direction}]"
