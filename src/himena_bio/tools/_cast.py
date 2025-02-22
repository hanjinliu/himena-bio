from __future__ import annotations

from typing import TYPE_CHECKING
from himena_bio.consts import SeqMeta

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


def cast_meta(meta) -> SeqMeta:
    if not isinstance(meta, SeqMeta):
        raise ValueError("Invalid metadata")
    return meta


def cast_seq_record(record) -> SeqRecord:
    from Bio.SeqRecord import SeqRecord

    if not isinstance(record, SeqRecord):
        raise ValueError("Invalid record")
    return record
