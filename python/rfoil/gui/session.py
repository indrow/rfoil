"""Application state: the analysis, accumulated polars and the worker thread.

Threading rule: while a job runs, the worker thread owns ``analysis``.  The
GUI thread only reads the :class:`Snapshot`, which is refreshed whenever the
session is idle and its inputs change.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PySide6.QtCore import QObject, QThread, Signal

from .. import Airfoil, Analysis, Paneling, Settings
from ..polar import Polar
from .worker import Job, Solution, Worker


def sequence(start: float, end: float, step: float) -> tuple[float, ...]:
    """Inclusive arithmetic sequence from ``start`` towards ``end``."""
    if step == 0 or start == end:
        return (start,)
    count = round(abs(end - start) / abs(step)) + 1
    return tuple(np.linspace(start, start + np.sign(end - start) * abs(step) * (count - 1), count))


@dataclass(frozen=True)
class Snapshot:
    """Analysis inputs and geometry, safe to read from the GUI thread."""

    airfoil: Airfoil
    settings: Settings
    paneling: Paneling
    panels: tuple[list[float], list[float]]
    shape: tuple[float, float, float, float]

    @classmethod
    def of(cls, analysis: Analysis) -> Snapshot:
        return cls(analysis.airfoil, analysis.settings, analysis.paneling, analysis.panels, analysis.shape())


class Session(QObject):
    """Owns the analysis; widgets talk to it and listen to its signals."""

    airfoil_changed = Signal()
    settings_changed = Signal()
    solved = Signal(object)  # Solution
    polars_changed = Signal()
    iteration = Signal(object)  # rfoil.Iteration
    progress = Signal(int, int)
    busy_changed = Signal(bool)
    message = Signal(str)

    _start = Signal(object, object)  # hands (analysis, job) to the worker thread

    def __init__(self) -> None:
        super().__init__()
        self.analysis = Analysis(Airfoil.naca("2412"), settings=Settings(reynolds=1e6))
        self.snapshot = Snapshot.of(self.analysis)
        self.polars: list[Polar] = []
        self.last: Solution | None = None
        self._busy = False
        self._target: Polar | None = None

        self._thread = QThread()
        self._thread.setObjectName("rfoil-solver")
        self._worker = Worker()
        self._worker.moveToThread(self._thread)
        self._start.connect(self._worker.run)
        self._worker.solved.connect(self._on_solved)
        self._worker.iteration.connect(self.iteration)
        self._worker.progress.connect(self.progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._thread.start()

    # ---------------------------------------------------------------- state

    @property
    def busy(self) -> bool:
        return self._busy

    @property
    def airfoil(self) -> Airfoil:
        return self.snapshot.airfoil

    @property
    def settings(self) -> Settings:
        return self.snapshot.settings

    @property
    def paneling(self) -> Paneling:
        return self.snapshot.paneling

    def panels(self) -> tuple[list[float], list[float]]:
        return self.snapshot.panels

    def shape(self) -> tuple[float, float, float, float]:
        return self.snapshot.shape

    def shutdown(self) -> None:
        self._worker.request_stop()
        self._thread.quit()
        self._thread.wait(5000)

    def _set_busy(self, busy: bool) -> None:
        self._busy = busy
        if not busy:
            self.snapshot = Snapshot.of(self.analysis)
        self.busy_changed.emit(busy)

    def _guard(self) -> bool:
        if self._busy:
            self.message.emit("Busy: stop the running analysis first")
        return not self._busy

    def _changed_geometry(self, text: str) -> None:
        self.snapshot = Snapshot.of(self.analysis)
        self.last = None
        self.airfoil_changed.emit()
        self.message.emit(text)

    # ------------------------------------------------------------- geometry

    def set_airfoil(self, airfoil: Airfoil) -> None:
        if self._guard():
            self.analysis.airfoil = airfoil
            self._changed_geometry(f"Loaded {airfoil.name} ({len(airfoil)} points)")

    def load_naca(self, designation: str) -> None:
        self.set_airfoil(Airfoil.naca(designation))

    def load_file(self, path: str | Path) -> None:
        path = Path(path)
        airfoil = Airfoil.parse(path.read_text(errors="replace"))
        if airfoil.name == "Unnamed":
            airfoil = airfoil.renamed(path.stem)
        self.set_airfoil(airfoil)

    def save_airfoil(self, path: str | Path) -> None:
        Path(path).write_text(self.airfoil.to_dat())
        self.message.emit(f"Saved coordinates to {path}")

    def set_paneling(self, paneling: Paneling) -> None:
        if self._guard():
            self.analysis.paneling = paneling
            self._changed_geometry(f"Repaneled with {paneling.nodes} nodes")

    # ------------------------------------------------------------- settings

    def set_settings(self, settings: Settings) -> None:
        if self._guard():
            self.analysis.settings = settings
            self.snapshot = Snapshot.of(self.analysis)
            self.settings_changed.emit()

    def reset_boundary_layer(self) -> None:
        if self._guard():
            self.analysis.reset_boundary_layer()
            self.message.emit("Boundary layer will be re-initialised on the next point")

    # ------------------------------------------------------------- analysis

    def run_alpha(self, start: float, end: float | None = None, step: float = 1.0) -> None:
        values = (start,) if end is None else sequence(start, end, step)
        self._run(Job("alpha", values))

    def run_cl(self, start: float, end: float | None = None, step: float = 0.1) -> None:
        values = (start,) if end is None else sequence(start, end, step)
        self._run(Job("cl", values))

    def stop(self) -> None:
        self._worker.request_stop()

    def _run(self, job: Job) -> None:
        if not self._guard():
            return
        self._target = self._polar_for(self.airfoil.name, self.settings)
        self._set_busy(True)
        kind = "α" if job.target == "alpha" else "CL"
        self.message.emit(f"Solving {len(job.values)} point(s) at prescribed {kind}…")
        self._start.emit(self.analysis, job)

    def _on_solved(self, solution: Solution) -> None:
        self.last = solution
        point = solution.point
        if point.converged and self._target is not None:
            self._target.add(point)
            self.polars_changed.emit()
        elif not point.converged:
            self.message.emit(f"α = {point.alpha:.2f}°: not converged (rms {point.rms_change:.1e})")
        self.solved.emit(solution)

    def _on_finished(self, text: str) -> None:
        self._set_busy(False)
        self.message.emit(text)

    def _on_failed(self, text: str) -> None:
        self._set_busy(False)
        self.message.emit(f"Error: {text}")

    # --------------------------------------------------------------- polars

    def _polar_for(self, airfoil: str, settings: Settings) -> Polar:
        for polar in self.polars:
            if polar.matches(airfoil, settings):
                return polar
        polar = Polar.for_settings(airfoil, settings)
        self.polars.append(polar)
        return polar

    def remove_polar(self, index: int) -> None:
        del self.polars[index]
        self.polars_changed.emit()

    def clear_polars(self) -> None:
        self.polars.clear()
        self.polars_changed.emit()

    def save_polar(self, index: int, path: str | Path) -> None:
        self.polars[index].save(path)
        self.message.emit(f"Saved polar to {path}")
