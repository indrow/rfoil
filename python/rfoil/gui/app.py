"""Entry point of the graphical application."""

from __future__ import annotations

import sys

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow, apply_initial_theme
from .session import Session


def main(argv: list[str] | None = None) -> int:
    app = QApplication(sys.argv if argv is None else argv)
    app.setApplicationName("rfoil")
    app.setOrganizationName("rfoil")
    apply_initial_theme(app)

    session = Session()
    window = MainWindow(session)
    for path in app.arguments()[1:]:
        session.load_file(path)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
