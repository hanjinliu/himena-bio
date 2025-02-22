from __future__ import annotations

from himena import WidgetDataModel
from qtpy import QtWidgets as QtW, QtGui
from himena.plugins import validate_protocol
from himena.consts import MonospaceFontFamily
from himena_bio.consts import Type
from Bio.Align import PairwiseAlignments


class QAlignmentView(QtW.QWidget):
    def __init__(self):
        super().__init__()
        self._ith = QtW.QSpinBox()
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
        self._ith.setMaximum(len(model.value) - 1)
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
        self._score.setText(f"Score = {self._alignments[index].score:.2f}")
        self._view.setPlainText(str(self._alignments[index]))
