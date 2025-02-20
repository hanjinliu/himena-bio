from __future__ import annotations

import re
from typing import TYPE_CHECKING, Callable
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from cmap import Color
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as _Tm

from himena import WidgetDataModel
from himena.consts import MonospaceFontFamily
from himena.plugins import validate_protocol
from himena_plasmid_editor.consts import Keys, ApeAnnotation

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
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self.setWordWrapMode(QtGui.QTextOption.WrapMode.WrapAnywhere)
        self.setLineWrapMode(QtW.QTextEdit.LineWrapMode.WidgetWidth)
        self.setUndoRedoEnabled(True)
        self.setAcceptRichText(False)
        self.setTabChangesFocus(True)
        self.setAcceptRichText(False)
        self.setAcceptDrops(True)
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
                    color = _parse_color(colors[0])
                    format.setBackground(color)
                    format.setForeground(
                        Qt.GlobalColor.white 
                        if color.lightness() < 128
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
            letter_annotations=self._seq_template.letter_annotations
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
    
    def _copy_selection(self, rev_comp: bool = False):
        if rev_comp:
            seq = self.textCursor().selectedText()
        else:
            seq = str(Seq(self.textCursor().selectedText()).reverse_complement())
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
        self._seq_edit = QSeqEdit(self)
        layout.addWidget(self._seq_edit)
        self._comment = QtW.QPlainTextEdit(self)
        self._comment.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self._comment.setWordWrapMode(QtGui.QTextOption.WrapMode.WordWrap)
        self._comment.setFixedHeight(80)
        layout.addWidget(self._comment)
        
        self._tm_method: Callable[[Seq], float] = _Tm.Tm_GC
        
        self._control = QSeqControl(self)
        self._seq_edit.selectionChanged.connect(self._selection_changed)

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        recs = model.value
        if len(recs) == 0:
            return
        choices: list[str] = []
        for rec in recs:
            if not isinstance(rec, SeqRecord):
                raise ValueError(f"Expected SeqRecord, got {type(rec)}")
            choices.append(rec.name)
        self._seq_choices.addItems(choices)
        self._set_record(recs[0])
    
    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        return WidgetDataModel(value=[self._seq_edit.to_record()])

    @validate_protocol
    def control_widget(self) -> QSeqControl:
        return self._control
    
    @validate_protocol
    def size_hint(self) -> tuple[int, int]:
        return 400, 300

    def _set_record(self, record: SeqRecord):
        self._seq_edit.set_record(record)
        cursor = self._seq_edit.textCursor()
        cursor.setPosition(0)
        self._seq_edit.setTextCursor(cursor)
        self._comment.setPlainText(record.annotations.get(ApeAnnotation.COMMENT, ""))

    def _selection_changed(self):
        cursor = self._seq_edit.textCursor()
        selection = cursor.selectedText()
        if selection:
            self._control._sel.set_value(f"{cursor.selectionStart()} - {cursor.selectionEnd()}")
            self._control._length.set_value(str(len(selection)))
            gc_count = sum(1 for nucleotide in selection if nucleotide in "GC")
            self._control._percent_gc.set_value(f"{(gc_count / len(selection)) * 100:.2f}%")
            tm = self._tm_method(Seq(selection))
            if tm < 0:
                self._control._tm.set_value(f"-- °C")
            else:
                self._control._tm.set_value(f"{tm:.1f}°C")
            _visible = True
        else:
            self._control._sel.set_value(f"{cursor.position() - 1} | {cursor.position()}")
            _visible = False
        self._control._length.set_visible(_visible)
        self._control._tm.set_visible(_visible)
        self._control._percent_gc.set_visible(_visible)

        
class QSeqControl(QtW.QWidget):
    def __init__(self, edit: QMultiSeqEdit):
        super().__init__()
        self._edit = edit
        layout = QtW.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        spacer = QtW.QWidget()
        layout.addWidget(spacer, stretch=10)
        self._sel = QVParam("Selection", 90)
        layout.addWidget(self._sel)
        self._length = QVParam("Length", 65)
        layout.addWidget(self._length)
        self._tm = QVParam("Tm", 65)
        layout.addWidget(self._tm)
        self._percent_gc = QVParam("%GC", 65)
        layout.addWidget(self._percent_gc)
        
        self._topology = QtW.QComboBox()
        self._topology.addItems(["linear", "circular"])
        self._topology.setFixedWidth(96)
        layout.addWidget(self._topology)

        
class QVParam(QtW.QWidget):
    def __init__(self, label: str, width: int = 96):
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
    
    def set_value(self, value: str):
        self._value.setText(value)

    def set_visible(self, visible: bool):
        self._label.setVisible(visible)
        self._value.setVisible(visible)

_GRAY_PATTERN = re.compile(r"gray(\d+)")

def _parse_color(color: str) -> QtGui.QColor:
    if match := _GRAY_PATTERN.match(color):
        val = round(255 * int(match.group(1)) / 100)
        return QtGui.QColor(val, val, val)
    return QtGui.QColor(Color(color).hex)
