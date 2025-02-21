from himena.plugins import register_widget_class
from himena_bio.consts import Type


def make_widget():
    from himena_bio.widgets.editor import QMultiSeqEdit

    return QMultiSeqEdit()


register_widget_class(Type.SEQS, make_widget)
