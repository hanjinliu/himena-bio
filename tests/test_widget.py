from himena.testing import choose_one_dialog_response
from himena import WidgetDataModel
from qtpy.QtCore import Qt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from himena_bio.widgets.editor import QSeqEdit, QMultiSeqEdit
from himena_bio.consts import Keys, Type
from pytestqt.qtbot import QtBot
import pytest

@pytest.fixture
def seq_control(qtbot: QtBot):
    """Test QSeqControl widget."""
    widget = QSeqEdit()
    qtbot.addWidget(widget)
    widget.set_keys_allowed(Keys.DNA)
    return widget

@pytest.fixture
def mseq_edit(qtbot: QtBot):
    widget = QMultiSeqEdit()
    qtbot.addWidget(widget)
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

def test_find_feature(himena_ui, seq_control: QSeqEdit):
    from Bio.SeqFeature import FeatureLocation

    record = SeqRecord(Seq("ACGTACGT"), id="test")
    feature1 = SeqFeature(FeatureLocation(0, 4), type="gene", id="feat1")
    feature2 = SeqFeature(FeatureLocation(4, 8), type="CDS", id="feat2")
    record.features = [feature1, feature2]
    seq_control.set_record(record)

    himena_ui.add_widget(seq_control)

    with choose_one_dialog_response(himena_ui, feature1):
        seq_control._find_feature()

def _selection(edit: QSeqEdit) -> tuple[int, int]:
    cursor = edit.textCursor()
    return cursor.selectionStart(), cursor.selectionEnd()

def test_finder_incremental(mseq_edit: QMultiSeqEdit):
    mseq_edit.update_model(
        WidgetDataModel(value=[SeqRecord(Seq("ttACGTaaACGAcgt"))], type=Type.DNA)
    )
    finder = mseq_edit._finder_widget
    edit = mseq_edit._seq_edit
    finder.show()
    finder.set_text("a")  # case insensitive
    assert _selection(edit) == (2, 3)
    finder.set_text("ac")  # should stay at the same position
    assert _selection(edit) == (2, 4)
    finder.set_text("acga")  # jump to the next match
    assert _selection(edit) == (8, 12)
    finder.set_text("ac")  # stay
    assert _selection(edit) == (8, 10)
    finder._find_next()
    assert _selection(edit) == (11, 13)
    finder._find_next()  # wrap around
    assert _selection(edit) == (2, 4)
    finder._find_prev()  # wrap around
    assert _selection(edit) == (11, 13)
    finder._find_prev()
    assert _selection(edit) == (8, 10)
    finder.set_text("ggg")  # not found
    assert _selection(edit) == (8, 10)

def test_finder_revcomp(mseq_edit: QMultiSeqEdit):
    mseq_edit.update_model(
        WidgetDataModel(value=[SeqRecord(Seq("AAAAGGCTTTTGCC"))], type=Type.DNA)
    )
    finder = mseq_edit._finder_widget
    edit = mseq_edit._seq_edit
    finder._also_find_reverse_complement.setChecked(True)
    finder.set_text("ggc")
    assert _selection(edit) == (4, 7)
    finder._find_next()
    assert _selection(edit) == (11, 14)  # GCC is the reverse complement

def _record_with_features() -> SeqRecord:
    rec = SeqRecord(Seq("AAAAACCCCCGGGGGTTTTT"), id="test")
    rec.features = [
        SeqFeature(SimpleLocation(0, 5, 1), type="f0"),
        SeqFeature(SimpleLocation(5, 10, -1), type="f1"),
        SeqFeature(SimpleLocation(12, 18, 1), type="f2"),
    ]
    return rec

def _feature_locs(edit: QSeqEdit):
    return [
        (f.type, int(f.location.start), int(f.location.end), f.location.strand)
        for f in edit._record.features
    ]

def test_delete_and_undo_features(seq_control: QSeqEdit):
    seq_control.set_record(_record_with_features())
    before = _feature_locs(seq_control)
    cursor = seq_control.textCursor()
    cursor.setPosition(4)
    cursor.setPosition(14, cursor.MoveMode.KeepAnchor)
    seq_control.delete_text(cursor)
    assert seq_control.toPlainText() == "AAAAGTTTTT"
    # f1 is completely deleted, f0 and f2 are clipped
    assert _feature_locs(seq_control) == [("f0", 0, 4, 1), ("f2", 4, 8, 1)]
    seq_control.undo()
    assert seq_control.toPlainText() == "AAAAACCCCCGGGGGTTTTT"
    assert _feature_locs(seq_control) == before
    seq_control.redo()
    assert _feature_locs(seq_control) == [("f0", 0, 4, 1), ("f2", 4, 8, 1)]

def test_replace_selection_and_undo(seq_control: QSeqEdit):
    seq_control.set_record(_record_with_features())
    before = _feature_locs(seq_control)
    cursor = seq_control.textCursor()
    cursor.setPosition(6)
    cursor.setPosition(8, cursor.MoveMode.KeepAnchor)
    seq_control.insert_text("TTTT", cursor)
    assert seq_control.toPlainText() == "AAAAACTTTTCCGGGGGTTTTT"
    assert _feature_locs(seq_control) == [
        ("f0", 0, 5, 1), ("f1", 5, 12, -1), ("f2", 14, 20, 1)
    ]
    seq_control.undo()
    assert seq_control.toPlainText() == "AAAAACCCCCGGGGGTTTTT"
    assert _feature_locs(seq_control) == before

def test_feature_delete_move_undo(seq_control: QSeqEdit):
    seq_control.set_record(_record_with_features())
    features = list(seq_control._record.features)

    def _types():
        return [f.type for f in seq_control._record.features]

    seq_control._delete_feature(features[1])
    assert _types() == ["f0", "f2"]
    seq_control.undo()
    assert _types() == ["f0", "f1", "f2"]
    seq_control._move_feature_front(features[0])
    assert _types() == ["f1", "f2", "f0"]
    seq_control.undo()
    assert _types() == ["f0", "f1", "f2"]
    seq_control._move_feature_back(features[2])
    assert _types() == ["f2", "f0", "f1"]
    seq_control.undo()
    assert _types() == ["f0", "f1", "f2"]

def test_backspace_delete_at_ends(seq_control: QSeqEdit, qtbot: QtBot):
    seq_control.setPlainText("ACGT")
    cursor = seq_control.textCursor()
    cursor.setPosition(0)
    seq_control.setTextCursor(cursor)
    qtbot.keyClick(seq_control, Qt.Key.Key_Backspace)
    assert seq_control.toPlainText() == "ACGT"
    cursor.setPosition(4)
    seq_control.setTextCursor(cursor)
    qtbot.keyClick(seq_control, Qt.Key.Key_Delete)
    assert seq_control.toPlainText() == "ACGT"

def test_input_not_modified(mseq_edit: QMultiSeqEdit):
    rec = _record_with_features()
    mseq_edit.update_model(WidgetDataModel(value=[rec], type=Type.DNA))
    assert not mseq_edit.is_modified()
    mseq_edit._seq_edit._delete_feature(mseq_edit._seq_edit._record.features[0])
    assert len(rec.features) == 3
    assert mseq_edit.is_modified()

def test_multi_records(mseq_edit: QMultiSeqEdit):
    recs = [
        SeqRecord(Seq("acgtacgt"), id="r0", name="r0"),
        SeqRecord(Seq("ggggcccc"), id="r1", name="r1"),
    ]
    recs[1].annotations["comment"] = "old comment\nApEinfo:methylated:1"
    mseq_edit.update_model(WidgetDataModel(value=recs, type=Type.SEQS))
    assert mseq_edit.model_type() == Type.DNA
    assert not mseq_edit.is_modified()
    assert mseq_edit._seq_choices.count() == 2
    mseq_edit._seq_edit.insert_text("T", mseq_edit._seq_edit.textCursor())
    mseq_edit._seq_choices.setCurrentIndex(1)
    assert mseq_edit._seq_edit.toPlainText() == "ggggcccc"
    assert mseq_edit._comment.toPlainText() == "old comment"
    mseq_edit._comment.setPlainText("new comment")
    model = mseq_edit.to_model()
    assert [str(r.seq) for r in model.value] == ["Tacgtacgt", "ggggcccc"]
    assert model.value[1].annotations["comment"] == (
        "new comment\nApEinfo:methylated:1"
    )
    assert model.metadata.current_index == 1
    # undo history is kept for each record
    mseq_edit._seq_choices.setCurrentIndex(0)
    mseq_edit._seq_edit.undo()
    assert mseq_edit._seq_edit.toPlainText() == "acgtacgt"
    # update_model again should not duplicate the entries
    mseq_edit.update_model(mseq_edit.to_model())
    assert mseq_edit._seq_choices.count() == 2

def test_protein_alignment(qtbot: QtBot):
    from himena_bio.tools.align import _pairwise_impl
    from himena_bio.widgets.alignment import QAlignmentView
    from Bio.Align import substitution_matrices

    seq0 = WidgetDataModel(value=[SeqRecord(Seq("mkvlaagiw"))], type=Type.PROTEIN)
    seq1 = WidgetDataModel(value=[SeqRecord(Seq("MKVLGIW"))], type=Type.PROTEIN)
    out = _pairwise_impl(
        seq0, seq1, mode="global",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10, extend_gap_score=-0.5,
    )
    view = QAlignmentView()
    qtbot.addWidget(view)
    view.update_model(out)
    assert "Score" in view._score.text()
    assert "MKVL" in view._view.toPlainText()
    for _ in range(5):
        view._ith._on_next()  # should not raise even if out of range
    assert view._ith.value() < len(out.value)
