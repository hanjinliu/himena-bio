from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Callable, Iterable
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.SeqUtils import MeltingTemp

from magicgui.widgets import Dialog
from cmap import Color
from himena import WidgetDataModel
from himena.widgets import set_status_tip, current_instance
from himena.types import is_subtype
from himena.qt.magicgui import get_type_map
from himena.consts import MonospaceFontFamily
from himena.plugins import validate_protocol
from himena.utils.collections import UndoRedoStack

from himena_bio.consts import Keys, ApeAnnotation, SeqMeta, Type
from himena_bio._utils import (
    feature_color,
    feature_to_slice,
    get_feature_label,
    map_location,
    copy_feature,
)
from himena_bio.widgets._feature_view import QFeatureView
from himena_bio.widgets._editor_finder import QSeqFinder
from himena_bio.widgets._base import char_to_qt_key, infer_seq_type
from himena_bio.widgets._editor_control import QSeqControl
from himena_bio.widgets import _editor_actions as _ea

_KEYS_MOVE = frozenset(
    [Qt.Key.Key_Left, Qt.Key.Key_Right, Qt.Key.Key_Up, Qt.Key.Key_Down,
    Qt.Key.Key_Home, Qt.Key.Key_End, Qt.Key.Key_PageUp, Qt.Key.Key_PageDown]
)  # fmt: skip

_MOLECULE_TYPES = {Type.DNA: "DNA", Type.RNA: "RNA", Type.PROTEIN: "protein"}


class QSeqEdit(QtW.QPlainTextEdit):
    hovered = QtCore.Signal(object)  # int or None
    edited = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self.setWordWrapMode(QtGui.QTextOption.WrapMode.WrapAnywhere)
        self.setLineWrapMode(QtW.QPlainTextEdit.LineWrapMode.WidgetWidth)
        # Undo/redo is implemented by the widget itself, not by the Qt document.
        self.setUndoRedoEnabled(False)
        self.setTabChangesFocus(True)
        self.setAcceptDrops(True)
        self.setMouseTracking(True)
        self._record = SeqRecord(Seq(""))
        self.set_keys_allowed(Keys.DNA, Keys.DNA_AMBIGUOUS)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        self._actions: dict[str, QtW.QAction] = {}
        self._setup_actions()
        self._undo_redo_stack = UndoRedoStack[_ea.EditorAction](size=100)

    def _setup_actions(self):
        """Setup actions for keyboard shortcuts and context menu."""

        def add_action(name: str, shortcut: str, callback):
            action = QtW.QAction(name, self)
            action.setShortcut(shortcut)
            action.triggered.connect(callback)
            action.setShortcutContext(Qt.ShortcutContext.WidgetShortcut)
            self._actions[name] = action
            self.addAction(action)

        add_action("Copy", "Ctrl+C", lambda: self._copy_selection(False))
        add_action("Cut", "Ctrl+X", lambda: self._cut_selection(False))
        add_action("Paste", "Ctrl+V", lambda: self._paste(False))
        add_action("Undo", "Ctrl+Z", self.undo)
        add_action("Redo", "Ctrl+Y", self.redo)
        add_action("Find", "Ctrl+F", self._open_finder)
        add_action("Copy (RevComp)", "Ctrl+Alt+C", lambda: self._copy_selection(True))
        add_action("Cut (RevComp)", "Ctrl+Alt+X", lambda: self._cut_selection(True))
        add_action("Paste (RevComp)", "Ctrl+Alt+V", lambda: self._paste(True))
        add_action("New Feature", "Alt+N", self._new_feature)
        add_action(
            "Edit Feature",
            "Alt+E",
            lambda: self._edit_feature(self._get_front_feature()),
        )
        add_action(
            "Delete Feature",
            "Alt+D",
            lambda: self._delete_feature(self._get_front_feature()),
        )
        add_action(
            "Move Feature Front",
            "Alt+F",
            lambda: self._move_feature_front(self._get_front_feature()),
        )
        add_action(
            "Move Feature Back",
            "Alt+B",
            lambda: self._move_feature_back(self._get_front_feature()),
        )
        add_action(
            "Find Feature",
            "Ctrl+Shift+F",
            self._find_feature,
        )

    def set_record(self, record: SeqRecord):
        """Set a new record. This resets the undo/redo history."""
        # copy the record so that editing does not affect the input.
        self._record = SeqRecord(
            record.seq,
            id=record.id,
            name=record.name,
            description=record.description,
            dbxrefs=list(record.dbxrefs),
            features=list(record.features),
            annotations=dict(record.annotations),
            letter_annotations=dict(record.letter_annotations),
        )
        self._undo_redo_stack = UndoRedoStack[_ea.EditorAction](size=100)
        self.setPlainText(str(record.seq))
        self.update_highlight()

    def update_highlight(self):
        """Update the background colors of the text to show the features."""
        doc_length = self.document().characterCount() - 1
        cursor = QtGui.QTextCursor(self.document())
        cursor.select(QtGui.QTextCursor.SelectionType.Document)
        cursor.setCharFormat(QtGui.QTextCharFormat())
        for feature in self._record.features:
            if feature.location is None:
                continue
            if (color := feature_color(feature)) is None:
                continue
            format = QtGui.QTextCharFormat()
            format.setBackground(color)
            r, g, b, _ = color.getRgbF()
            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            format.setForeground(
                Qt.GlobalColor.white if luminance < 0.5 else Qt.GlobalColor.black
            )
            for part in feature.location.parts:
                start = min(max(int(part.start), 0), doc_length)
                end = min(max(int(part.end), 0), doc_length)
                if start >= end:
                    continue
                cursor.setPosition(start)
                cursor.setPosition(end, QtGui.QTextCursor.MoveMode.KeepAnchor)
                cursor.setCharFormat(format)

    def to_record(self) -> SeqRecord:
        seq = Seq(self.toPlainText())
        letter_annotations = {
            k: v
            for k, v in self._record.letter_annotations.items()
            if len(v) == len(seq)  # dropped if the sequence length changed
        }
        return SeqRecord(
            seq=seq,
            id=self._record.id,
            name=self._record.name,
            description=self._record.description,
            dbxrefs=list(self._record.dbxrefs),
            features=list(self._record.features),
            annotations=dict(self._record.annotations),
            letter_annotations=letter_annotations,
        )

    def keyPressEvent(self, event: QtGui.QKeyEvent):
        _mod = event.modifiers()
        _key = event.key()
        if _key in _KEYS_MOVE:
            return super().keyPressEvent(event)
        if _mod & Qt.KeyboardModifier.ControlModifier:
            _alt_on = bool(_mod & Qt.KeyboardModifier.AltModifier)
            if _key == Qt.Key.Key_C:
                self._copy_selection(_alt_on)
            elif _key == Qt.Key.Key_X:
                self._cut_selection(_alt_on)
            elif _key == Qt.Key.Key_A:
                self.selectAll()
            elif _key == Qt.Key.Key_V:
                self._paste(_alt_on)
            elif _key == Qt.Key.Key_Z and not _alt_on:
                self.undo()
            elif _key == Qt.Key.Key_Y and not _alt_on:
                self.redo()
            return
        elif _mod & Qt.KeyboardModifier.AltModifier:
            return None
        else:
            if _key in self._keys_allowed:
                _key_char = event.text()
                self.insert_text(_key_char, self.textCursor())
            elif _key == Qt.Key.Key_Backspace:
                self._backspace_event()
            elif _key == Qt.Key.Key_Delete:
                self._delete_event()
            return None

    def mouseMoveEvent(self, e):
        self.hovered.emit(self.cursorForPosition(e.pos()).position())
        return super().mouseMoveEvent(e)

    def leaveEvent(self, a0):
        self.hovered.emit(None)
        return super().leaveEvent(a0)

    def _copy_selection(self, rev_comp: bool = False):
        seq = self.textCursor().selectedText()
        if seq == "":
            return
        if rev_comp:
            seq = str(Seq(seq).reverse_complement())
        clipboard = QtW.QApplication.clipboard()
        clipboard.setText(seq)

    def _paste(self, rev_comp: bool = False):
        clipboard = QtW.QApplication.clipboard()
        # ignore spaces, line breaks and numbers (copied from GenBank-like text)
        seq = "".join(c for c in clipboard.text() if not (c.isspace() or c.isdigit()))
        if seq == "":
            return
        if invalid := {c for c in seq if c.upper() not in self._chars_allowed}:
            _set_status_tip(f"Cannot paste invalid characters: {sorted(invalid)}")
            return
        if rev_comp:
            seq = str(Seq(seq).reverse_complement())
        self.insert_text(seq, self.textCursor())

    def _cut_selection(self, rev_comp: bool = False):
        cursor = self.textCursor()
        if not cursor.hasSelection():
            return
        self._copy_selection(rev_comp)
        self.delete_text(cursor)

    def _after_edit(self):
        self.update_highlight()
        self.edited.emit()

    def insert_text(
        self,
        text: str,
        cursor: QtGui.QTextCursor,
        record_undo: bool = True,
    ):
        """Insert the text using the given text cursor."""
        if text == "" and not cursor.hasSelection():
            return
        actions: list[_ea.EditorAction] = []
        if cursor.hasSelection():
            actions.extend(self._delete_selection(cursor))
        start = cursor.position()
        nchars = len(text)
        for index, feature in enumerate(self._record.features):
            feature_shifted = _feature_after_insert(feature, start, nchars)
            if feature_shifted is not feature:
                action = _ea.EditFeatureAction(
                    index=index, old=feature, new=feature_shifted
                )
                action.apply(self)
                actions.append(action)
        actions.append(_ea.InsertSeqAction(pos=start, seq=text))
        cursor.insertText(text)
        self._after_edit()
        if record_undo:
            self._undo_redo_stack.push(_ea.CompositeAction(actions))

    def delete_text(self, cursor: QtGui.QTextCursor, record_undo: bool = True):
        """Delete the selected text using the given text cursor."""
        if not cursor.hasSelection():
            return
        actions = self._delete_selection(cursor)
        self._after_edit()
        if record_undo:
            self._undo_redo_stack.push(_ea.CompositeAction(actions))

    def _delete_selection(self, cursor: QtGui.QTextCursor) -> list[_ea.EditorAction]:
        """Delete the selected text and update features, return the actions."""
        start, end = cursor.selectionStart(), cursor.selectionEnd()
        size = self.document().characterCount() - 1
        actions: list[_ea.EditorAction] = [
            _ea.DeleteSeqAction(
                start=start, length=end - start, seq=cursor.selectedText()
            )
        ]
        # iterate in the reverse order so that deleting a feature does not affect
        # the indices of the remaining ones.
        for index in reversed(range(len(self._record.features))):
            feature = self._record.features[index]
            feature_new = _feature_after_delete(feature, start, end, size)
            if feature_new is not feature:
                action = _ea.EditFeatureAction(
                    index=index, old=feature, new=feature_new
                )
                action.apply(self)
                actions.append(action)
        cursor.removeSelectedText()
        return actions

    def undo(self):
        """Undo the last edit action."""
        if action := self._undo_redo_stack.undo():
            action.invert().apply(self)
            self._after_edit()
        return None

    def redo(self):
        """Redo the last undone edit action."""
        if action := self._undo_redo_stack.redo():
            action.apply(self)
            self._after_edit()
        return None

    def _open_finder(self):
        finder: QSeqFinder = self.parentWidget()._finder_widget
        if text := self.textCursor().selectedText():
            finder.set_text(text)
        finder.show()

    def _backspace_event(self):
        cursor = self.textCursor()
        if not cursor.hasSelection():
            if cursor.position() == 0:
                return
            cursor.movePosition(
                QtGui.QTextCursor.MoveOperation.PreviousCharacter,
                QtGui.QTextCursor.MoveMode.KeepAnchor,
            )
        self.delete_text(cursor)

    def _delete_event(self):
        cursor = self.textCursor()
        if not cursor.hasSelection():
            if cursor.atEnd():
                return
            cursor.movePosition(
                QtGui.QTextCursor.MoveOperation.NextCharacter,
                QtGui.QTextCursor.MoveMode.KeepAnchor,
            )
        self.delete_text(cursor)

    def _make_context_menu(self) -> QtW.QMenu:
        menu = QtW.QMenu()
        menu.addAction(self._actions["Cut"])
        menu.addAction(self._actions["Cut (RevComp)"])
        menu.addAction(self._actions["Copy"])
        menu.addAction(self._actions["Copy (RevComp)"])
        menu.addAction(self._actions["Paste"])
        menu.addAction(self._actions["Paste (RevComp)"])
        menu.addSeparator()
        menu.addAction(self._actions["New Feature"])
        menu.addAction(self._actions["Edit Feature"])
        menu.addAction(self._actions["Delete Feature"])
        menu.addAction(self._actions["Move Feature Front"])
        menu.addAction(self._actions["Move Feature Back"])

        can_paste = QtW.QApplication.clipboard().text() != ""
        self._actions["Paste"].setEnabled(can_paste)
        self._actions["Paste (RevComp)"].setEnabled(can_paste)
        return menu

    def _show_context_menu(self, pos: QtCore.QPoint):
        cursor_pos = self.cursorForPosition(pos).position()

        cursor = self.textCursor()
        start, end = cursor.selectionStart(), cursor.selectionEnd()
        if not start <= cursor_pos < end:
            cursor.setPosition(cursor_pos)
            self.setTextCursor(cursor)
        menu = self._make_context_menu()
        menu.exec(self.mapToGlobal(pos))

    def _features_under_pos(self, pos: int) -> list[SeqFeature]:
        return [
            feat
            for feat in self._record.features
            if feat.location is not None and pos in feat
        ]

    def _new_feature(self):
        cursor = self.textCursor()
        start, end = cursor.selectionStart(), cursor.selectionEnd()
        if start == end:
            return
        kwargs = self._feature_qualifiers_from_dialog()
        if kwargs is None:
            return
        feature = SeqFeature(location=SimpleLocation(start, end, strand=1), **kwargs)
        action = _ea.EditFeatureAction(
            index=len(self._record.features), old=None, new=feature
        )
        action.apply(self)
        self._undo_redo_stack.push(action)
        self._after_edit()

    def _get_front_feature(self) -> SeqFeature | None:
        """Return the feature under the cursor that is displayed at the front."""
        pos = self.textCursor().position()
        if features := self._features_under_pos(pos):
            return features[-1]
        return None

    def _edit_feature(self, feature: SeqFeature | None):
        if feature is None:
            return
        kwargs = self._feature_qualifiers_from_dialog(
            name=get_feature_label(feature),
            type=feature.type,
            fcolor=feature.qualifiers.get(ApeAnnotation.FWCOLOR, ["cyan"])[0],
            rcolor=feature.qualifiers.get(ApeAnnotation.RVCOLOR, ["cyan"])[0],
        )
        if kwargs is None:
            return
        index = self._record.features.index(feature)
        feature_new = copy_feature(feature)
        feature_new.type = kwargs["type"]
        feature_new.qualifiers.update(kwargs["qualifiers"])
        action = _ea.EditFeatureAction(index=index, old=feature, new=feature_new)
        action.apply(self)
        self._undo_redo_stack.push(action)
        self._after_edit()

    def _delete_feature(self, feature: SeqFeature | None):
        """Delete the feature from the current record."""
        if feature is None:
            return
        index = self._record.features.index(feature)
        action = _ea.EditFeatureAction(index=index, old=feature, new=None)
        action.apply(self)
        self._undo_redo_stack.push(action)
        self._after_edit()

    def _move_feature_front(self, feature: SeqFeature | None):
        """Make sure the feature is visible by moving it the last of the list."""
        if feature is None:
            return
        self._move_feature(feature, len(self._record.features) - 1)

    def _move_feature_back(self, feature: SeqFeature | None):
        if feature is None:
            return
        self._move_feature(feature, 0)

    def _move_feature(self, feature: SeqFeature, new: int):
        old = self._record.features.index(feature)
        if old == new:
            return
        action = _ea.MoveFeatureAction(old=old, new=new)
        action.apply(self)
        self._undo_redo_stack.push(action)
        self._after_edit()

    def _find_feature(self):
        if not self._record.features:
            _set_status_tip("No feature in this sequence.")
            return
        choices: list[tuple[str, SeqFeature]] = []
        for feat in self._record.features:
            txt = f"{get_feature_label(feat)} (type: {feat.type}, location: {feat.location})"
            choices.append((txt, feat))
        if feat := current_instance().exec_choose_one_dialog(
            message="Choose a feature",
            choices=choices,
            how="palette",
        ):
            self._select_feature(feat, 0)

    def _select_feature(self, feature: SeqFeature, nth: int):
        try:
            start, end = feature_to_slice(feature, nth)
        except NotImplementedError:
            return
        cursor = self.textCursor()
        cursor.setPosition(max(start, 0))
        cursor.setPosition(
            min(end, self.document().characterCount() - 1),
            QtGui.QTextCursor.MoveMode.KeepAnchor,
        )
        self.setTextCursor(cursor)

    def _feature_qualifiers_from_dialog(
        self,
        name: str = "Unnamed",
        type: str = "misc_feature",
        fcolor: str = "cyan",
        rcolor: str = "cyan",
    ) -> dict[str, Any] | None:
        typemap = get_type_map()
        w_name = typemap.create_widget(value=name, label="Name")
        w_type = typemap.create_widget(value=type, label="Type")
        w_fcolor = typemap.create_widget(value=Color(fcolor), label="Forward Color")
        w_rcolor = typemap.create_widget(value=Color(rcolor), label="Reverse Color")
        dlg = Dialog(widgets=[w_name, w_type, w_fcolor, w_rcolor])
        if dlg.exec():
            return {
                "type": w_type.value or "misc_feature",
                "qualifiers": {
                    ApeAnnotation.LABEL: [w_name.value],
                    ApeAnnotation.FWCOLOR: [Color(w_fcolor.value).hex],
                    ApeAnnotation.RVCOLOR: [Color(w_rcolor.value).hex],
                },
            }
        return None

    def set_keys_allowed(
        self,
        keys: Iterable[str],
        paste_keys: Iterable[str] | None = None,
    ):
        """Set characters allowed for typing (`keys`) and pasting (`paste_keys`)."""
        keys = frozenset(keys)
        self._keys_allowed = frozenset(char_to_qt_key(char) for char in keys)
        self._chars_allowed = keys if paste_keys is None else frozenset(paste_keys)


class QMultiSeqEdit(QtW.QWidget):
    """Sequence editor widget."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtW.QVBoxLayout(self)

        self._seq_choices = QtW.QComboBox(self)
        self._feature_view = QFeatureView(self)
        self._seq_edit = QSeqEdit(self)
        self._finder_widget = QSeqFinder(self, self._seq_edit)

        layout.addWidget(self._seq_choices)
        self._seq_choices.hide()
        self._seq_choices.currentIndexChanged.connect(self._on_choice_changed)

        self._feature_view.setFixedHeight(40)
        self._feature_view.clicked.connect(self._on_view_clicked)
        self._feature_view.hovered.connect(self._on_view_hovered)
        layout.addWidget(self._feature_view)

        layout.addWidget(self._finder_widget)
        self._finder_widget.hide()

        self._seq_edit.hovered.connect(self._seq_edit_hovered)
        self._seq_edit.edited.connect(self._seq_edited)
        layout.addWidget(self._seq_edit)

        self._comment = QtW.QPlainTextEdit(self)
        self._comment.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self._comment.setWordWrapMode(QtGui.QTextOption.WrapMode.WordWrap)
        self._comment.setFixedHeight(80)
        self._comment.setPlaceholderText("Comment")
        self._comment.textChanged.connect(self._set_modified)
        layout.addWidget(self._comment)

        self._tm_method: Callable[[Seq], float] = MeltingTemp.Tm_GC

        self._control = QSeqControl(self)
        self._control._topology.changed.connect(self._set_modified)
        self._seq_edit.selectionChanged.connect(self._selection_changed)
        self._seq_edit.cursorPositionChanged.connect(self._selection_changed)
        self._model_type = Type.DNA
        self._extension_default = ".ape"

        self._records: list[SeqRecord] = []
        self._undo_stacks: list[UndoRedoStack[_ea.EditorAction]] = []
        self._current_index = -1
        self._modified = False
        self._is_updating = False

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        recs = model.value
        if isinstance(recs, (str, Seq, SeqRecord)):
            recs = [recs]
        recs_normed = [_to_seq_record(rec) for rec in recs]
        if len(recs_normed) == 0:
            recs_normed = [SeqRecord(Seq(""))]
        if model.type == Type.SEQS:
            _mtype = infer_seq_type(str(recs_normed[0].seq))
        else:
            _mtype = model.type
        if is_subtype(_mtype, Type.DNA):
            self._seq_edit.set_keys_allowed(Keys.DNA, Keys.DNA_AMBIGUOUS)
        elif is_subtype(_mtype, Type.RNA):
            self._seq_edit.set_keys_allowed(Keys.RNA, Keys.RNA_AMBIGUOUS)
        elif is_subtype(_mtype, Type.PROTEIN):
            self._seq_edit.set_keys_allowed(Keys.PROTEIN)
        else:
            raise NotImplementedError(f"Unsupported model type: {_mtype}")
        self._model_type = _mtype
        if ext := model.extension_default:
            self._extension_default = ext

        self._records = recs_normed
        self._undo_stacks = [UndoRedoStack(size=100) for _ in recs_normed]
        index = 0
        if isinstance(meta := model.metadata, SeqMeta):
            if 0 <= meta.current_index < len(recs_normed):
                index = meta.current_index
        with self._updating():
            self._seq_choices.clear()
            self._seq_choices.addItems([_record_name(rec) for rec in recs_normed])
            self._seq_choices.setCurrentIndex(index)
        self._seq_choices.setVisible(len(recs_normed) > 1)
        self._switch_to(index)
        if isinstance(meta, SeqMeta):
            self._set_selection(*meta.selection)
        self._modified = False

    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        cursor = self._seq_edit.textCursor()
        self._store_current()
        return WidgetDataModel(
            value=list(self._records),
            type=self._model_type,
            metadata=SeqMeta(
                current_index=max(self._current_index, 0),
                selection=(cursor.selectionStart(), cursor.selectionEnd()),
            ),
            extension_default=self._extension_default,
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

    @validate_protocol
    def is_modified(self) -> bool:
        return self._modified

    @validate_protocol
    def set_modified(self, value: bool) -> None:
        self._modified = value

    @validate_protocol
    def widget_added_callback(self):
        self._feature_view.auto_range()

    def setFocus(self):
        self._seq_edit.setFocus()

    @contextmanager
    def _updating(self):
        """Programmatic update (should not be considered as user modification)."""
        was_updating = self._is_updating
        self._is_updating = True
        try:
            yield
        finally:
            self._is_updating = was_updating

    def _set_modified(self, *_):
        if not self._is_updating:
            self._modified = True

    def _current_record(self) -> SeqRecord:
        """Build the record from the current state of the widgets."""
        record = self._seq_edit.to_record()
        record.annotations["topology"] = self._control._topology.value
        if mol_type := _MOLECULE_TYPES.get(self._model_type):
            record.annotations.setdefault("molecule_type", mol_type)
        old_comment = record.annotations.get(ApeAnnotation.COMMENT, "")
        if isinstance(old_comment, list):
            old_comment = "\n".join(old_comment)
        comment = _restore_ape_meta(self._comment.toPlainText(), old_comment)
        if comment:
            record.annotations[ApeAnnotation.COMMENT] = comment
        else:
            record.annotations.pop(ApeAnnotation.COMMENT, None)
        return record

    def _store_current(self):
        if 0 <= self._current_index < len(self._records):
            self._records[self._current_index] = self._current_record()

    def _on_choice_changed(self, index: int):
        if self._is_updating or index < 0:
            return
        self._store_current()
        self._switch_to(index)

    def _switch_to(self, index: int):
        self._current_index = index
        self._set_record(self._records[index])
        # each record has its own undo/redo history
        self._seq_edit._undo_redo_stack = self._undo_stacks[index]

    def _set_record(self, record: SeqRecord):
        with self._updating():
            self._seq_edit.set_record(record)
            self._feature_view.set_record(record)
            self._feature_view.auto_range()
            self._set_selection(0, 0)
            comment = record.annotations.get(ApeAnnotation.COMMENT, "")
            if isinstance(comment, list):
                comment = "\n".join(comment)
            self._comment.setPlainText(_remove_ape_meta(comment))
            self._control._topology.set_value(
                record.annotations.get("topology", "linear")
            )

    def _set_selection(self, start: int, end: int):
        size = self._seq_edit.document().characterCount() - 1
        cursor = self._seq_edit.textCursor()
        cursor.setPosition(min(max(start, 0), size))
        cursor.setPosition(
            min(max(end, 0), size), QtGui.QTextCursor.MoveMode.KeepAnchor
        )
        self._seq_edit.setTextCursor(cursor)

    def _is_nucleotide(self) -> bool:
        return is_subtype(self.model_type(), Type.DNA) or is_subtype(
            self.model_type(), Type.RNA
        )

    def _selection_changed(self):
        cursor = self._seq_edit.textCursor()
        selection = cursor.selectedText()
        offset = 1 if self._control._is_one_start.isChecked() else 0
        is_nucleotide = self._is_nucleotide()
        if selection:
            # 0-start: Python slice [start, end), 1-start: ApE style [start, end]
            self._control._sel.set_value(
                f"{cursor.selectionStart() + offset} - {cursor.selectionEnd()}"
            )
            self._control._length.set_value(str(len(selection)))
            # nucleotide specific info
            if is_nucleotide:
                selection_upper = selection.upper()
                gc_count = sum(1 for nuc in selection_upper if nuc in "GCS")
                self._control._percent_gc.set_value(
                    f"{(gc_count / len(selection)) * 100:.2f}%"
                )
                try:
                    tm = self._tm_method(Seq(selection_upper))
                except Exception:
                    tm = -1
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
        self._control._tm.set_visible(_visible and is_nucleotide)
        self._control._percent_gc.set_visible(_visible and is_nucleotide)

        if is_nucleotide and len(selection) >= 3:
            # translate the selection and show in status tip
            adjusted_length = len(selection) // 3 * 3
            try:
                amino_acids = str(Seq(selection[:adjusted_length]).translate())
            except Exception:
                amino_acids = ""
            if len(amino_acids) > 20:
                amino_acids = amino_acids[:10] + "..." + amino_acids[-10:]
            _set_status_tip(amino_acids, duration=5)

        has_selection = cursor.hasSelection()
        self._seq_edit._actions["Cut"].setEnabled(has_selection)
        self._seq_edit._actions["Cut (RevComp)"].setEnabled(has_selection)
        self._seq_edit._actions["Copy"].setEnabled(has_selection)
        self._seq_edit._actions["Copy (RevComp)"].setEnabled(has_selection)
        self._seq_edit._actions["New Feature"].setEnabled(has_selection)

        has_feature = len(self._seq_edit._features_under_pos(cursor.position())) > 0
        for name in [
            "Edit Feature",
            "Delete Feature",
            "Move Feature Front",
            "Move Feature Back",
        ]:
            self._seq_edit._actions[name].setEnabled(has_feature)

    def _seq_edit_hovered(self, pos: int | None):
        if pos is None:
            self._control._feature_label.setText("")
            self._control._hover_pos.set_value("")
            return
        offset = 1 if self._control._is_one_start.isChecked() else 0
        self._control._hover_pos.set_value(str(pos + offset))
        feature_labels = [
            get_feature_label(feat) for feat in self._seq_edit._features_under_pos(pos)
        ]
        tip = ", ".join(feature_labels)
        self._control._feature_label.setText(tip)

    def _seq_edited(self):
        self._set_modified()
        self._feature_view.set_record(self._seq_edit.to_record())
        self._selection_changed()

    def _on_view_clicked(self, feature: SeqFeature | None, nth: int):
        if feature is None:
            self._set_selection(nth, nth)
            return
        self._seq_edit._select_feature(feature, nth)

    def _on_view_hovered(self, feature: SeqFeature | None, pos: int | None):
        if isinstance(feature, SeqFeature):
            label = get_feature_label(feature)
            self._control._feature_label.setText(label)
        else:
            self._control._feature_label.setText("")
        if pos is None:
            self._control._hover_pos.set_value("")
        else:
            offset = 1 if self._control._is_one_start.isChecked() else 0
            self._control._hover_pos.set_value(str(pos + offset))


def _to_seq_record(rec) -> SeqRecord:
    if isinstance(rec, SeqRecord):
        return rec
    return SeqRecord(Seq(str(rec)))


def _record_name(rec: SeqRecord) -> str:
    for name in [rec.name, rec.id]:
        if name and not name.startswith("<unknown"):
            return name
    return "Untitled"


def _set_status_tip(msg: str, duration: float = 5):
    try:
        set_status_tip(msg, duration=duration)
    except Exception:  # no main window
        pass


def _feature_after_insert(feature: SeqFeature, pos: int, nchars: int) -> SeqFeature:
    """Return the feature updated for text insertion (same object if unchanged)."""
    if (loc := feature.location) is None or nchars == 0:
        return feature
    changed = False
    new_parts: list[SimpleLocation] = []
    for part in loc.parts:
        start, end = int(part.start), int(part.end)
        if start >= pos:
            start += nchars
            changed = True
        if end > pos:
            end += nchars
            changed = True
        new_parts.append(SimpleLocation(start, end, strand=part.strand))
    if not changed:
        return feature
    if len(new_parts) == 1:
        return copy_feature(feature, new_parts[0])
    return copy_feature(
        feature, CompoundLocation(new_parts, operator=getattr(loc, "operator", "join"))
    )


def _feature_after_delete(
    feature: SeqFeature, start: int, end: int, size: int
) -> SeqFeature | None:
    """Return the feature updated for text deletion.

    Returns the same object if unchanged, and None if the feature is completely
    deleted.
    """
    if (loc := feature.location) is None:
        return feature
    if all(int(part.end) <= start for part in loc.parts):
        return feature
    loc_new = map_location(loc, [(0, start), (end, size)])
    if loc_new is None:
        return None
    return copy_feature(feature, loc_new)


def _remove_ape_meta(comment: str) -> str:
    lines = comment.splitlines()
    if lines and lines[-1].startswith("ApEinfo:methylated:"):
        lines.pop(-1)
    return "\n".join(lines)


def _restore_ape_meta(comment: str, old_comment: str) -> str:
    """Add the ApE meta line removed by `_remove_ape_meta` back to the comment."""
    old_lines = old_comment.splitlines()
    if old_lines and old_lines[-1].startswith("ApEinfo:methylated:"):
        if comment:
            return f"{comment}\n{old_lines[-1]}"
        return old_lines[-1]
    return comment
