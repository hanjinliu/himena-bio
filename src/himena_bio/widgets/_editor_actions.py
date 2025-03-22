from __future__ import annotations

from typing import TYPE_CHECKING
from dataclasses import dataclass
from qtpy import QtGui
from Bio.SeqFeature import SeqFeature

if TYPE_CHECKING:
    from typing import Self
    from himena_bio.widgets.editor import QSeqEdit


class EditorAction:
    def apply(self, widget: QSeqEdit):
        raise NotImplementedError

    def invert(self) -> Self:
        raise NotImplementedError


@dataclass
class EditFeatureAction(EditorAction):
    """Edit a feature, including adding/deleting one."""

    index: int
    old: SeqFeature | None
    new: SeqFeature | None

    def apply(self, widget: QSeqEdit):
        if self.new is None:
            widget._record.features.pop(self.index)
        elif len(widget._record.features) <= self.index:
            widget._record.features.append(self.new)
        else:
            widget._record.features[self.index] = self.new

    def invert(self) -> Self:
        return EditFeatureAction(index=self.index, old=self.new, new=self.old)


@dataclass
class MoveFeatureAction(EditorAction):
    old: int
    new: int

    def apply(self, widget: QSeqEdit):
        if self.old < self.new:
            new = self.new - 1
        else:
            new = self.new
        widget._record.features.insert(new, widget._record.features.pop(self.old))

    def invert(self) -> Self:
        return MoveFeatureAction(old=self.new, new=self.old)


@dataclass
class InsertSeqAction(EditorAction):
    pos: int
    seq: str

    def apply(self, widget: QSeqEdit):
        cursor = widget.textCursor()
        cursor.setPosition(self.pos)
        widget.insert_text(self.seq, cursor, record_undo=False)

    def invert(self) -> Self:
        return DeleteSeqAction(start=self.pos, length=len(self.seq), seq=self.seq)


@dataclass
class DeleteSeqAction(EditorAction):
    start: int
    length: int
    seq: str

    def apply(self, widget: QSeqEdit):
        cursor = widget.textCursor()
        cursor.setPosition(self.start)
        cursor.setPosition(
            self.start + self.length, QtGui.QTextCursor.MoveMode.KeepAnchor
        )
        widget.delete_text(cursor, record_undo=False)

    def invert(self) -> Self:
        return InsertSeqAction(pos=self.start, seq=self.seq)


@dataclass
class CompositeAction(EditorAction):
    actions: list[EditorAction]

    def apply(self, widget: QSeqEdit):
        for action in self.actions:
            action.apply(widget)

    def invert(self) -> Self:
        return CompositeAction(
            actions=[action.invert() for action in reversed(self.actions)]
        )
