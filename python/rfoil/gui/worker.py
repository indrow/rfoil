"""Runs analyses off the GUI thread."""

from __future__ import annotations

import threading
from dataclasses import dataclass
from typing import Literal

from PySide6.QtCore import QObject, Signal, Slot

from .. import Analysis, BoundaryLayerSide, OperatingPoint, Surface

#: Consecutive non-converged points after which a sequence stops.
FAILURE_LIMIT = 4


@dataclass(frozen=True)
class Job:
    """A sequence of operating points with α (degrees) or CL prescribed."""

    target: Literal["alpha", "cl"]
    values: tuple[float, ...]


@dataclass(frozen=True)
class Solution:
    """Everything needed to plot one solved operating point."""

    point: OperatingPoint
    surface: Surface
    boundary_layer: tuple[BoundaryLayerSide, BoundaryLayerSide] | None
    wake: tuple[list[float], list[float]]


class Worker(QObject):
    """Lives in a worker thread and owns the analysis while a job runs."""

    solved = Signal(object)  # Solution
    iteration = Signal(object)  # rfoil.Iteration
    progress = Signal(int, int)  # done, total
    finished = Signal(str)
    failed = Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self._stop = threading.Event()

    def request_stop(self) -> None:
        self._stop.set()

    @Slot(object, object)
    def run(self, analysis: Analysis, job: Job) -> None:
        self._stop.clear()
        solve = analysis.solve_alpha if job.target == "alpha" else analysis.solve_cl
        failures = converged = 0
        try:
            for done, value in enumerate(job.values):
                if self._stop.is_set():
                    break
                point = solve(value, progress=self.iteration.emit)
                self.solved.emit(_snapshot(analysis, point))
                self.progress.emit(done + 1, len(job.values))
                if point.converged:
                    failures, converged = 0, converged + 1
                else:
                    failures += 1
                    if failures >= FAILURE_LIMIT:
                        self.finished.emit(f"Stopped after {failures} unconverged points")
                        return
        except (RuntimeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        verb = "Stopped" if self._stop.is_set() else "Done"
        self.finished.emit(f"{verb}: {converged} converged point(s)")


def _snapshot(analysis: Analysis, point: OperatingPoint) -> Solution:
    return Solution(
        point=point,
        surface=analysis.surface(),
        boundary_layer=analysis.boundary_layer() if point.viscous else None,
        wake=analysis.wake,
    )
