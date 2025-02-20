from __future__ import annotations


from typing import Callable
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from superqt import QElidingLabel
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as _Tm

from himena import WidgetDataModel
from himena.widgets import set_status_tip
from himena.qt.magicgui._toggle_switch import QLabeledToggleSwitch
from himena.consts import MonospaceFontFamily
from himena.plugins import validate_protocol
from himena_plasmid_editor.consts import Keys, ApeAnnotation, SeqMeta, Type
from himena_plasmid_editor._utils import parse_ape_color, get_feature_label
from himena_plasmid_editor.widgets._feature_view import QFeatureView


def _char_to_qt_key(char: str) -> Qt.Key:
    return getattr(Qt.Key, f"Key_{char.upper()}")


_KEYS_DNA = frozenset(_char_to_qt_key(char) for char in Keys.DNA)
_KEYS_RNA = frozenset(_char_to_qt_key(char) for char in Keys.RNA)
_KEYS_DEL = frozenset([Qt.Key.Key_Backspace, Qt.Key.Key_Delete])
_KEYS_MOVE = frozenset(
    [Qt.Key.Key_Left, Qt.Key.Key_Right, Qt.Key.Key_Up, Qt.Key.Key_Down,
    Qt.Key.Key_Home, Qt.Key.Key_End, Qt.Key.Key_PageUp, Qt.Key.Key_PageDown]
)  # fmt: skip
_FEATURE = QtGui.QTextFormat.Property.UserProperty
_LOC_ID = QtGui.QTextFormat.Property.UserProperty + 1


class _FeatureBox:
    def __init__(self, value: SeqFeature):
        self.value = value


class QSeqEdit(QtW.QTextEdit):
    hovered = QtCore.Signal(object)  # int or None

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self.setWordWrapMode(QtGui.QTextOption.WrapMode.WrapAnywhere)
        self.setLineWrapMode(QtW.QTextEdit.LineWrapMode.WidgetWidth)
        self.setUndoRedoEnabled(True)
        self.setTabChangesFocus(True)
        self.setAcceptRichText(False)
        self.setAcceptDrops(True)
        self.setMouseTracking(True)
        self._seq_template = SeqRecord(Seq(""))
        self._keys_allowed = _KEYS_DNA | _KEYS_DEL | _KEYS_MOVE

    def set_record(self, record: SeqRecord):
        self.setPlainText(str(record.seq))
        cursor = self.textCursor()
        format = QtGui.QTextCharFormat()
        for feature in record.features:
            if isinstance(loc := feature.location, SimpleLocation):
                parts = [(int(loc.start), int(loc.end))]
            elif isinstance(loc := feature.location, CompoundLocation):
                parts = [(int(sub.start), int(sub.end)) for sub in loc]
            else:
                parts = []
            for ith, (start, end) in enumerate(parts):
                cursor.setPosition(start)
                cursor.setPosition(end, QtGui.QTextCursor.MoveMode.KeepAnchor)
                if colors := feature.qualifiers.get(ApeAnnotation.FWCOLOR):
                    color = parse_ape_color(colors[0])
                    format.setBackground(color)
                    r, g, b, _ = color.getRgbF()
                    luminance = 0.299 * r + 0.587 * g + 0.114 * b
                    format.setForeground(
                        Qt.GlobalColor.white
                        if luminance < 0.5
                        else Qt.GlobalColor.black
                    )
                    format.setProperty(_FEATURE, _FeatureBox(feature))
                    format.setProperty(_LOC_ID, ith)
                    cursor.setCharFormat(format)
        self._seq_template = record
        return None

    def to_record(self) -> SeqRecord:
        # TODO: update features
        features = self._seq_template.features
        return SeqRecord(
            seq=self.toPlainText(),
            id=self._seq_template.id,
            name=self._seq_template.name,
            description=self._seq_template.description,
            dbxrefs=self._seq_template.dbxrefs,
            features=features,
            annotations=self._seq_template.annotations,
            letter_annotations=self._seq_template.letter_annotations,
        )

    def keyPressEvent(self, event: QtGui.QKeyEvent):
        _mod = event.modifiers()
        _key = event.key()
        if _mod & Qt.KeyboardModifier.ControlModifier:
            _shift_on = bool(_mod & Qt.KeyboardModifier.ShiftModifier)
            if _key == Qt.Key.Key_C:
                self._copy_selection(_shift_on)
            elif _key == Qt.Key.Key_X:
                self._cut_selection(_shift_on)
            elif _key == Qt.Key.Key_V:
                self._paste(_shift_on)
            return
        elif _mod & Qt.KeyboardModifier.AltModifier:
            return None
        else:
            if _key in self._keys_allowed:
                return super().keyPressEvent(event)
            return None

    def mouseMoveEvent(self, e):
        self.hovered.emit(self.cursorForPosition(e.pos()).position())
        return super().mouseMoveEvent(e)

    def leaveEvent(self, a0):
        self.hovered.emit(None)
        return super().leaveEvent(a0)

    def _copy_selection(self, rev_comp: bool = False):
        if rev_comp:
            seq = str(Seq(self.textCursor().selectedText()).reverse_complement())
        else:
            seq = self.textCursor().selectedText()
        clipboard = QtW.QApplication.clipboard()
        clipboard.setText(seq)

    def _paste(self, rev_comp: bool = False):
        clipboard = QtW.QApplication.clipboard()
        seq = clipboard.text()
        if rev_comp:
            seq = str(Seq(seq).reverse_complement())
        self.textCursor().insertText(seq)

    def _cut_selection(self, rev_comp: bool = False):
        self._copy_selection(rev_comp)
        self.textCursor().removeSelectedText()

    def _post_text_insert(self, start: int, text: str):
        nchars = len(text)
        cursor = self.textCursor()
        cursor.setPosition(start)
        fmt = cursor.charFormat()
        if isinstance(feature := fmt.property(_FEATURE), _FeatureBox):
            if isinstance(loc := feature.value.location, SimpleLocation):
                loc.start += nchars
                loc.end += nchars
            elif isinstance(loc := feature.value.location, CompoundLocation):
                loc_id = fmt.property(_LOC_ID)
                loc[loc_id].start += nchars
                loc[loc_id].end += nchars
            feature.location = loc

    def _post_text_remove(self, start: int, nchars: int):
        cursor = self.textCursor()
        cursor.setPosition(start)
        fmt = cursor.charFormat()
        if isinstance(feature := fmt.property(_FEATURE), _FeatureBox):
            if isinstance(loc := feature.value.location, SimpleLocation):
                loc.start -= nchars
                loc.end -= nchars
            elif isinstance(loc := feature.value.location, CompoundLocation):
                loc_id = fmt.property(_LOC_ID)
                loc[loc_id].start -= nchars
                loc[loc_id].end -= nchars
            feature.location = loc


class QMultiSeqEdit(QtW.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtW.QVBoxLayout(self)

        self._seq_choices = QtW.QComboBox(self)
        layout.addWidget(self._seq_choices)

        self._feature_view = QFeatureView()
        self._feature_view.setFixedHeight(40)
        self._feature_view.clicked.connect(self._on_view_clicked)
        self._feature_view.hovered.connect(self._on_view_hovered)
        layout.addWidget(self._feature_view)

        self._seq_edit = QSeqEdit(self)
        self._seq_edit.hovered.connect(self._seq_edit_hovered)
        layout.addWidget(self._seq_edit)

        self._comment = QtW.QPlainTextEdit(self)
        self._comment.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self._comment.setWordWrapMode(QtGui.QTextOption.WrapMode.WordWrap)
        self._comment.setFixedHeight(80)
        layout.addWidget(self._comment)

        self._tm_method: Callable[[Seq], float] = _Tm.Tm_GC

        self._control = QSeqControl(self)
        self._seq_edit.selectionChanged.connect(self._selection_changed)
        self._model_type = Type.DNA

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        recs = model.value
        if len(recs) == 0:
            return
        choices: list[str] = []
        for rec in recs:
            if not isinstance(rec, SeqRecord):
                rec = SeqRecord(Seq(rec))
            choices.append(rec.name)
        self._seq_choices.addItems(choices)
        self._set_record(recs[0])
        self._seq_choices.setVisible(len(recs) > 1)
        self._model_type = model.type
        if model.type == Type.DNA:
            self._seq_edit._keys_allowed = frozenset(
                _char_to_qt_key(char) for char in Keys.DNA
            )
        elif model.type == Type.RNA:
            self._seq_edit._keys_allowed = frozenset(
                _char_to_qt_key(char) for char in Keys.RNA
            )
        elif model.type == Type.PROTEIN:
            self._seq_edit._keys_allowed = frozenset(
                _char_to_qt_key(char) for char in Keys.PROTEIN
            )
        else:
            raise NotImplementedError(f"Unsupported model type: {model.type}")

    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        cursor = self._seq_edit.textCursor()
        return WidgetDataModel(
            value=[self._seq_edit.to_record()],
            type=self._model_type,
            metadata=SeqMeta(
                current_index=self._seq_choices.currentIndex(),
                selection=(cursor.selectionStart(), cursor.selectionEnd()),
            ),
        )

    @validate_protocol
    def model_type(self) -> str:
        return self._model_type

    @validate_protocol
    def control_widget(self) -> QSeqControl:
        return self._control

    @validate_protocol
    def size_hint(self) -> tuple[int, int]:
        return 400, 400

    def _set_record(self, record: SeqRecord):
        self._seq_edit.set_record(record)
        self._feature_view.set_record(record)
        cursor = self._seq_edit.textCursor()
        cursor.setPosition(0)
        self._seq_edit.setTextCursor(cursor)
        self._comment.setPlainText(record.annotations.get(ApeAnnotation.COMMENT, ""))

    def _selection_changed(self):
        cursor = self._seq_edit.textCursor()
        selection = cursor.selectedText()
        offset = 1 if self._control._is_one_start.isChecked() else 0
        if selection:
            self._control._sel.set_value(
                f"{cursor.selectionStart() + offset} - {cursor.selectionEnd() + offset}"
            )
            self._control._length.set_value(str(len(selection)))
            gc_count = sum(1 for nucleotide in selection if nucleotide in "GC")
            self._control._percent_gc.set_value(
                f"{(gc_count / len(selection)) * 100:.2f}%"
            )
            tm = self._tm_method(Seq(selection))
            if tm < 0:
                self._control._tm.set_value("-- °C")
            else:
                self._control._tm.set_value(f"{tm:.1f}°C")
            _visible = True
        else:
            pos = cursor.position() + offset
            self._control._sel.set_value(f"{pos - 1} | {pos}")
            _visible = False
        self._control._length.set_visible(_visible)
        self._control._tm.set_visible(_visible)
        self._control._percent_gc.set_visible(_visible)

        # translate the selection and show in status tip
        adjusted_length = len(selection) // 3 * 3
        amino_acids = str(Seq(selection[:adjusted_length]).translate())
        if len(amino_acids) > 20:
            amino_acids = amino_acids[:10] + "..." + amino_acids[-10:]
        set_status_tip(amino_acids, duration=5)

    def _seq_edit_hovered(self, pos: int | None):
        if pos is None:
            self._control._feature_label.setText("")
            self._control._hover_pos.set_value("")
            return
        offset = 1 if self._control._is_one_start.isChecked() else 0
        self._control._hover_pos.set_value(str(pos + offset))
        features = [
            get_feature_label(feat)
            for feat in self._seq_edit._seq_template.features
            if pos in feat
        ]
        tip = ", ".join(features)
        self._control._feature_label.setText(tip)

    def _on_view_clicked(self, feature: SeqFeature | None, nth: int):
        if feature is None:
            cursor = self._seq_edit.textCursor()
            cursor.setPosition(nth)
            self._seq_edit.ensureCursorVisible()
            self._seq_edit.setTextCursor(cursor)
            return
        if isinstance(loc := feature.location, SimpleLocation):
            start, end = int(loc.start), int(loc.end)
        elif isinstance(loc := feature.location, CompoundLocation):
            start, end = int(loc[nth].start), int(loc[nth].end)
        else:
            return
        cursor = self._seq_edit.textCursor()
        cursor.setPosition(start)
        self._seq_edit.ensureCursorVisible()
        cursor.setPosition(end, QtGui.QTextCursor.MoveMode.KeepAnchor)
        self._seq_edit.setTextCursor(cursor)

    def _on_view_hovered(self, feature: SeqFeature | None, nth: int):
        if isinstance(feature, SeqFeature):
            label = get_feature_label(feature)
            self._control._feature_label.setText(label)
        else:
            self._control._feature_label.setText("")


class QSeqControl(QtW.QWidget):
    def __init__(self, edit: QMultiSeqEdit):
        super().__init__()
        self._edit = edit
        layout = QtW.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._feature_label = QElidingLabel()
        self._feature_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self._feature_label, stretch=10)
        self._is_one_start = QLabeledToggleSwitch()
        self._is_one_start.setText("1-start")
        self._is_one_start.setToolTip(
            "Toggle between 0-based (Python style) and 1-based (ApE style) numbering for base positions."
        )
        self._is_one_start.toggled.connect(edit._selection_changed)
        layout.addWidget(self._is_one_start)

        self._hover_pos = QVParam("Pos", width=36, tooltip="Hovered position")
        layout.addWidget(self._hover_pos)

        self._sel = QVParam(
            "Selection", width=70, tooltip="Current text cursor selection"
        )
        layout.addWidget(self._sel)

        self._length = QVParam(
            "Length", width=45, tooltip="Length of the current text cursor selection"
        )
        layout.addWidget(self._length)

        self._tm = QVParam("Tm", width=40, tooltip="The calculated melting temperature")
        layout.addWidget(self._tm)
        self._percent_gc = QVParam(
            "%GC",
            width=45,
            tooltip="Percentage of GC content at the text cursor selection.",
        )
        layout.addWidget(self._percent_gc)

        self._topology = QtW.QComboBox()
        self._topology.setToolTip("Topology of the sequence.")
        self._topology.addItems(["linear", "circular"])
        self._topology.setFixedWidth(72)
        layout.addWidget(self._topology)


class QVParam(QtW.QWidget):
    def __init__(self, label: str, width: int = 96, tooltip: str = ""):
        super().__init__()
        layout = QtW.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(1)
        self._label = QtW.QLabel(label)
        self._value = QtW.QLabel()
        font = self._label.font()
        font.setPointSize(8)
        font.setBold(True)
        self._label.setFont(font)
        layout.addWidget(self._label)
        layout.addWidget(self._value)
        self.setFixedWidth(width)
        self.setToolTip(tooltip)

    def set_value(self, value: str):
        self._value.setText(value)

    def set_visible(self, visible: bool):
        self._label.setVisible(visible)
        self._value.setVisible(visible)
