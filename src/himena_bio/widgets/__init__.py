from himena.plugins import register_widget_class
from himena_bio.consts import Type


def multi_seq_edit():
    from himena_bio.widgets.editor import QMultiSeqEdit

    return QMultiSeqEdit()


multi_seq_edit.__himena_widget_id__ = "himena-bio:QMultiSeqEdit"
multi_seq_edit.__himena_display_name__ = "Multi-sequence editor"


def alignment_view():
    from himena_bio.widgets.alignment import QAlignmentView

    return QAlignmentView()


alignment_view.__himena_widget_id__ = "himena-bio:QAlignmentView"
alignment_view.__himena_display_name__ = "Alignment viewer"


def ab1_view():
    from himena_bio.widgets.ab1 import QAB1View

    return QAB1View()


ab1_view.__himena_widget_id__ = "himena-bio:QAB1View"
ab1_view.__himena_display_name__ = "AB1 viewer"

register_widget_class(Type.SEQS, multi_seq_edit)
register_widget_class(Type.ALIGNMENT, alignment_view)
register_widget_class(Type.DNA_ABI, ab1_view)
