from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING
from Bio import SeqIO, AlignIO

from himena import WidgetDataModel
from himena.plugins import register_reader_plugin, register_writer_plugin
from himena_bio.consts import Type

if TYPE_CHECKING:
    from Bio.SeqIO import SeqRecord


@register_reader_plugin
def read_fasta(path: Path):
    """Read a FASTA file."""
    return WidgetDataModel(value=_seq_parse(path, "fasta"), type=Type.SEQS)


@read_fasta.define_matcher
def _(path: Path):
    return _matcher_impl(path, [".fasta", ".fa", ".seq"], Type.SEQS)


@register_reader_plugin
def read_gb(path: Path):
    """Read a GenBank file."""
    return WidgetDataModel(value=_seq_parse(path, "genbank"), type=Type.DNA)


@read_gb.define_matcher
def _(path: Path):
    return _matcher_impl(path, [".gb", ".gbk", ".ape"], Type.DNA)


@register_writer_plugin
def write_dna(path: Path, model: WidgetDataModel):
    if path.suffix in [".gb", ".gbk", ".ape"] and model.type == Type.DNA:
        if len(model.value) == 1:
            SeqIO.write(model.value, path, "genbank")
        else:
            raise ValueError("Only one sequence can be written to a GenBank file.")
    elif path.suffix in [".fasta", ".fa"]:
        SeqIO.write(model.value, path, "fasta")
    else:
        raise ValueError("Unsupported file format.")


@write_dna.define_matcher
def _(path: Path, data: WidgetDataModel):
    return data.type == Type.SEQS


@register_reader_plugin
def read_alignment(path: Path):
    """Read a sequence alignment file."""
    if path.suffix in [".aln", ".clustal", ".clw", ".clustalw"]:
        fmt = "clustal"
    elif path.suffix in [".phy"]:
        fmt = "phylip"
    elif path.suffix in [".nex", ".nexus"]:
        fmt = "nexus"
    elif path.suffix in [".emboss"]:
        fmt = "emboss"
    elif path.suffix in [".stk"]:
        fmt = "stockholm"

    return WidgetDataModel(value=list(AlignIO.parse(path, fmt)), type=Type.ALIGNMENT)


@read_alignment.define_matcher
def _(path: Path):
    return _matcher_impl(
        path,
        [".aln", ".phy", ".clustal", ".clw", ".clustalw", ".nex", ".nexus", ".emboss",
         ".stk"],
        Type.ALIGNMENT,
    )  # fmt: skip


@register_reader_plugin
def read_ab1(path: Path):
    """Read an AB1 file."""
    return WidgetDataModel(value=_seq_parse(path, "abi"), type=Type.DNA_ABI)


@read_ab1.define_matcher
def _(path: Path):
    return _matcher_impl(path, [".ab1"], Type.DNA_ABI)


def _seq_parse(path: Path, format: str) -> list[SeqRecord]:
    return list(SeqIO.parse(path, format))


def _matcher_impl(path: Path, allowed: list[str], typ: str) -> str | None:
    if path.suffix in allowed:
        return typ
    return None
