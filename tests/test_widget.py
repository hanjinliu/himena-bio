from himena_bio.widgets.editor import QSeqEdit
from himena_bio.consts import Keys
from pytestqt.qtbot import QtBot
import pytest

@pytest.fixture
def seq_control(qtbot: QtBot):
    """Test QSeqControl widget."""
    widget = QSeqEdit()
    qtbot.addWidget(widget)
    widget.set_keys_allowed(Keys.DNA)
    return widget

def test_type_text(seq_control: QSeqEdit, qtbot: QtBot):
    seq_control.setPlainText("ACGT")
    assert seq_control.toPlainText() == "ACGT"
    cursor = seq_control.textCursor()
    cursor.setPosition(2)
    seq_control.setTextCursor(cursor)
    qtbot.keyClick(seq_control, "q")  # not allowed
    assert seq_control.toPlainText() == "ACGT"
    seq_control.undo()
    assert seq_control.toPlainText() == "ACGT"
    seq_control.redo()
    assert seq_control.toPlainText() == "ACGT"

    qtbot.keyClick(seq_control, "t")  # allowed
    assert seq_control.toPlainText() == "ACtGT"
    seq_control.undo()
    assert seq_control.toPlainText() == "ACGT"
    seq_control.redo()
    assert seq_control.toPlainText() == "ACtGT"
