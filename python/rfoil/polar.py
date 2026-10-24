"""Accumulated polars and their file formats."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._core import OperatingPoint, Settings

#: Columns of a polar, in the order of the plain-text polar file.
COLUMNS = ("alpha", "cl", "cd", "cdp", "cm", "xtr_top", "xtr_bottom")


@dataclass
class Polar:
    """Converged operating points for one airfoil and one set of conditions."""

    airfoil: str
    reynolds: float
    mach: float
    ncrit: float
    xtr_top: float = 1.0
    xtr_bottom: float = 1.0
    points: list[OperatingPoint] = field(default_factory=list)

    @classmethod
    def for_settings(cls, airfoil: str, settings: Settings) -> Polar:
        return cls(
            airfoil=airfoil,
            reynolds=settings.reynolds,
            mach=settings.mach,
            ncrit=settings.ncrit,
            xtr_top=settings.xtr_top,
            xtr_bottom=settings.xtr_bottom,
        )

    @property
    def label(self) -> str:
        if self.reynolds <= 0:
            flow = "inviscid"
        else:
            flow = f"Re {self.reynolds / 1e6:g}M"
        if self.mach > 0:
            flow += f"  M {self.mach:g}"
        return f"{self.airfoil}  {flow}"

    def matches(self, airfoil: str, settings: Settings) -> bool:
        """True if points computed with ``settings`` belong in this polar."""
        return (
            self.airfoil == airfoil
            and self.reynolds == settings.reynolds
            and self.mach == settings.mach
            and self.ncrit == settings.ncrit
            and self.xtr_top == settings.xtr_top
            and self.xtr_bottom == settings.xtr_bottom
        )

    def add(self, point: OperatingPoint) -> None:
        """Inserts a converged point, replacing any existing one at the same α."""
        if not point.converged:
            return
        self.points = [p for p in self.points if abs(p.alpha - point.alpha) > 1e-6]
        self.points.append(point)
        self.points.sort(key=lambda p: p.alpha)

    def column(self, name: str) -> np.ndarray:
        return np.array([getattr(p, name) for p in self.points], dtype=float)

    def __len__(self) -> int:
        return len(self.points)

    def to_text(self) -> str:
        """The polar as a plain-text table with a header block."""
        lines = [
            "",
            "       rfoil polar",
            "",
            f" Calculated polar for: {self.airfoil}",
            "",
            f" xtrf =   {self.xtr_top:.3f} (top)        {self.xtr_bottom:.3f} (bottom)",
            f" Mach = {self.mach:7.3f}     Re = {self.reynolds / 1e6:9.3f} e 6     Ncrit = {self.ncrit:7.3f}",
            "",
            "   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr",
            "  ------ -------- --------- --------- -------- -------- --------",
        ]
        for p in self.points:
            lines.append(
                f"  {p.alpha:6.3f} {p.cl:8.4f} {p.cd:9.5f} {p.cdp:9.5f}"
                f" {p.cm:8.4f} {p.xtr_top:8.4f} {p.xtr_bottom:8.4f}"
            )
        return "\n".join(lines) + "\n"

    def to_csv(self) -> str:
        rows = [",".join(COLUMNS)]
        rows += [",".join(f"{getattr(p, c):.6g}" for c in COLUMNS) for p in self.points]
        return "\n".join(rows) + "\n"

    def save(self, path: str | Path) -> None:
        path = Path(path)
        text = self.to_csv() if path.suffix.lower() == ".csv" else self.to_text()
        path.write_text(text)
