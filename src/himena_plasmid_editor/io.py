from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterator
from Bio import SeqIO

from himena import WidgetDataModel
from himena.plugins import register_reader_plugin
from himena_plasmid_editor.consts import Type

if TYPE_CHECKING:
    from Bio.SeqIO import SeqRecord

@register_reader_plugin
def read_fasta(path: Path):
    return WidgetDataModel(value=list(_seq_parse(path, "fasta")), type=Type.SEQS)

@read_fasta.define_matcher
def _(path: Path):
    return _matcher_impl(path, (".fasta", ".fa"))

@register_reader_plugin
def read_gb(path: Path):
    return WidgetDataModel(value=list(_seq_parse(path, "genbank")), type=Type.SEQS)

@read_gb.define_matcher
def _(path: Path):
    return _matcher_impl(path, (".gb", ".gbk", ".ape"))

def _seq_parse(path, format: str) -> Iterator[SeqRecord]:
    yield from SeqIO.parse(str(path), format)

def _matcher_impl(path: Path, allowed) -> str | None:
    if path.suffix in allowed:
        return Type.SEQS
    return None
