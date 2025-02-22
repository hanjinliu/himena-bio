from __future__ import annotations

from qtpy import QtWidgets as QtW
from qtpy.QtCore import Qt
from himena_bio.consts import Type


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


def char_to_qt_key(char: str) -> Qt.Key:
    if char.isalnum():
        return getattr(Qt.Key, f"Key_{char.upper()}")
    if char == " ":
        return Qt.Key.Key_Space
    if char == "*":
        return Qt.Key.Key_Asterisk
    if char == "-":
        return Qt.Key.Key_Minus
    if char == "@":
        return Qt.Key.Key_At
    raise NotImplementedError(f"Unsupported character: {char}")


def infer_seq_type(seq: str) -> str:
    if set(seq) <= set("ATGCN"):
        return Type.DNA
    elif set(seq) <= set("AUGCNYRWSKMBDHV"):
        return Type.RNA
    elif set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"):
        return Type.PROTEIN
    else:
        raise ValueError("Unsupported sequence type.")
