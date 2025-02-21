from __future__ import annotations

from qtpy import QtWidgets as QtW
from qtpy.QtCore import Qt


class QBaseGraphicsScene(QtW.QGraphicsScene):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._grab_source: QtW.QGraphicsItem | None = None

    def grabSource(self) -> QtW.QGraphicsItem | None:
        return self._grab_source

    def setGrabSource(self, item: QtW.QGraphicsItem | None):
        self._grab_source = item


class QBaseGraphicsView(QtW.QGraphicsView):
    def __init__(self):
        scene = QBaseGraphicsScene()
        super().__init__(scene)
        self.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignHCenter)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setMouseTracking(True)

    def addItem(self, item: QtW.QGraphicsItem):
        self.scene().addItem(item)

    def scene(self) -> QBaseGraphicsScene:
        return super().scene()
