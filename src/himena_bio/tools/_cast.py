from __future__ import annotations

from typing import TYPE_CHECKING
from himena_bio.consts import SeqMeta

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord
    from himena import WidgetDataModel


def cast_meta(meta) -> SeqMeta:
    if not isinstance(meta, SeqMeta):
        raise ValueError("Invalid metadata")
    return meta


def cast_seq_record(record) -> SeqRecord:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if isinstance(record, (str, Seq)):
        return SeqRecord(Seq(str(record)))
    if not isinstance(record, SeqRecord):
        raise ValueError("Invalid record")
    return record


def current_record(model: WidgetDataModel) -> SeqRecord:
    """Get the record currently shown in the widget (the first one if unknown)."""
    records = model.value
    if len(records) == 0:
        raise ValueError(f"{model.title!r} has no sequence.")
    index = 0
    if isinstance(meta := model.metadata, SeqMeta):
        if 0 <= meta.current_index < len(records):
            index = meta.current_index
    return cast_seq_record(records[index])
