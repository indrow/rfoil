"""Airfoil selection, coordinate files and paneling."""

from __future__ import annotations

from PySide6.QtWidgets import QFileDialog, QLineEdit, QPushButton, QWidget

from ... import Paneling
from ..session import Session
from ..widgets import Card, ValueLabel, hbox, int_spin, spin

FILE_FILTER = "Airfoil coordinates (*.dat *.txt *.xy);;All files (*)"


class AirfoilPanel(Card):
    def __init__(self, session: Session, parent: QWidget | None = None) -> None:
        super().__init__("Airfoil", parent)
        self.session = session

        self.naca = QLineEdit(placeholderText="NACA 4/5-digit, e.g. 4412")
        generate = QPushButton("Generate")
        generate.clicked.connect(self._generate)
        self.naca.returnPressed.connect(self._generate)
        self.form.addRow(hbox(self.naca, generate))

        open_btn = QPushButton("Open…")
        open_btn.clicked.connect(self.open_file)
        save_btn = QPushButton("Save…")
        save_btn.clicked.connect(self.save_file)
        self.form.addRow(hbox(open_btn, save_btn))

        self.name = ValueLabel("{}")
        self.thickness = ValueLabel("{}")
        self.camber = ValueLabel("{}")
        self.row("Name", self.name)
        self.row("Thickness", self.thickness)
        self.row("Camber", self.camber)

        p = session.paneling
        self.nodes = self.row("Panel nodes", int_spin(p.nodes, 40, 400, 10), "Number of panel nodes")
        self.bunching = self.row(
            "LE bunching", spin(p.bunching, 0.0, 3.0, 0.1), "Curvature attraction: 0 = uniform"
        )
        self.te_ratio = self.row(
            "TE/LE density", spin(p.te_le_ratio, 0.01, 2.0, 0.05), "TE / LE panel density ratio"
        )
        repanel = QPushButton("Repanel")
        repanel.clicked.connect(self._repanel)
        self.form.addRow(repanel)

        self._editors = (self.naca, generate, open_btn, repanel, self.nodes, self.bunching, self.te_ratio)
        session.airfoil_changed.connect(self.refresh)
        session.busy_changed.connect(lambda busy: [w.setEnabled(not busy) for w in self._editors])
        self.refresh()

    def refresh(self) -> None:
        thickness, x_thickness, camber, x_camber = self.session.shape()
        self.name.show_value(self.session.airfoil.name)
        self.thickness.show_value(f"{100 * thickness:.2f} %  at x = {x_thickness:.3f}")
        self.camber.show_value(f"{100 * camber:.2f} %  at x = {x_camber:.3f}")

    def _generate(self) -> None:
        text = self.naca.text().strip()
        if text:
            self._attempt(lambda: self.session.load_naca(text))

    def open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open airfoil", "", FILE_FILTER)
        if path:
            self._attempt(lambda: self.session.load_file(path))

    def save_file(self) -> None:
        name = self.session.airfoil.name.replace(" ", "_") + ".dat"
        path, _ = QFileDialog.getSaveFileName(self, "Save airfoil", name, FILE_FILTER)
        if path:
            self._attempt(lambda: self.session.save_airfoil(path))

    def _repanel(self) -> None:
        paneling = Paneling(
            nodes=self.nodes.value(),
            bunching=self.bunching.value(),
            te_le_ratio=self.te_ratio.value(),
        )
        self._attempt(lambda: self.session.set_paneling(paneling))

    def _attempt(self, action) -> None:
        try:
            action()
        except (ValueError, OSError) as exc:
            self.session.message.emit(f"Error: {exc}")
