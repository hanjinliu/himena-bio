from __future__ import annotations

from himena import WidgetDataModel
from qtpy import QtWidgets as QtW, QtGui, QtCore
from himena.plugins import validate_protocol
from himena.consts import MonospaceFontFamily
from himena_bio.consts import Type
from Bio.Align import PairwiseAlignments


class QAlignmentView(QtW.QWidget):
    def __init__(self):
        super().__init__()
        self._ith = QAlignmentSpinBox()
        self._score = QtW.QLabel()
        self._view = QtW.QPlainTextEdit()
        self._view.setReadOnly(True)
        self._view.setWordWrapMode(QtGui.QTextOption.WrapMode.NoWrap)
        self._view.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        layout = QtW.QVBoxLayout(self)
        layout.addWidget(self._ith)
        layout.addWidget(self._view)
        self._ith.valueChanged.connect(self._on_index_changed)
        self._alignments: PairwiseAlignments | None = None
        self._model_type = Type.ALIGNMENT

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        if not isinstance(model.value, PairwiseAlignments):
            raise ValueError("Invalid alignment type")
        self._alignments = model.value
        self._model_type = model.type
        self._ith.setValue(0)
        self._on_index_changed(0)

    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        return WidgetDataModel(value=self._alignments, type=self.model_type())

    @validate_protocol
    def model_type(self) -> str:
        return self._model_type

    @validate_protocol
    def size_hint(self) -> tuple[int, int]:
        return 420, 500

    def _on_index_changed(self, index: int):
        try:
            aln = self._alignments[index]
        except StopIteration:
            self._ith.setMaximum(index)
            return
        self._score.setText(f"Score = {aln.score:.2f}")
        self._view.setPlainText(str(aln))


class QAlignmentSpinBox(QtW.QWidget):
    valueChanged = QtCore.Signal(int)

    def __init__(self):
        super().__init__()
        layout = QtW.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._left = QtW.QPushButton("◀")
        self._left.clicked.connect(self._on_prev)
        self._left.setFixedWidth(30)
        self._left.setEnabled(False)
        self._label = QtW.QLabel("0")
        self._label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self._right = QtW.QPushButton("▶")
        self._right.clicked.connect(self._on_next)
        self._right.setFixedWidth(30)

        layout.addWidget(self._left)
        layout.addWidget(self._label)
        layout.addWidget(self._right)

        self._max_value = float("inf")

    def _on_next(self):
        if self.value() < self._max_value:
            self._label.setText(str(int(self._label.text()) + 1))
            if self.value() > 0:
                self._left.setEnabled(True)
            self.valueChanged.emit(self.value())

    def _on_prev(self):
        if self.value() > 0:
            self._label.setText(str(int(self._label.text()) - 1))
            if self.value() == 0:
                self._left.setEnabled(False)
            if self.value() < self._max_value:
                self._right.setEnabled(True)
            self.valueChanged.emit(self.value())

    def value(self) -> int:
        return int(self._label.text())

    def setValue(self, value: int):
        self._label.setText(str(value))
        self._left.setEnabled(value > 0)
        self._right.setEnabled(value < self._max_value)

    def setMaximum(self, value: int):
        self._max_value = value
        if self.value() > value:
            self._label.setText(str(value))
        self._right.setEnabled(self.value() < value)
