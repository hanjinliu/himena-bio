from types import SimpleNamespace
from himena.standards.model_meta import BaseMetadata


class Type(SimpleNamespace):
    SEQS = "bio-seqs"
    DNA = "bio-seqs.dna"
    RNA = "bio-seqs.rna"
    PROTEIN = "bio-seqs.protein"


class Keys(SimpleNamespace):
    DNA = frozenset(["A", "T", "C", "G", "N"])
    DNA_AMBIGUOUS = frozenset(
        ["A", "T", "C", "G", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    )
    RNA = frozenset(["A", "U", "C", "G", "N"])
    RNA_AMBIGUOUS = frozenset(
        ["A", "U", "C", "G", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    )
    PROTEIN = frozenset(
        [
            "A",
            "R",
            "N",
            "D",
            "C",
            "Q",
            "E",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
            "B",
            "Z",
            "X",
            "*",
        ]
    )


class ApeAnnotation(SimpleNamespace):
    FWCOLOR = "ApEinfo_fwdcolor"
    RVCOLOR = "ApEinfo_revcolor"
    COMMENT = "comment"


class SeqMeta(BaseMetadata):
    current_index: int
    selection: tuple[int, int]
