from typing import Literal
from himena import Parametric, WidgetDataModel
from himena.plugins import register_function, configure_gui
from himena_bio.consts import Type
from himena_bio.tools._cast import current_record

_PROTEIN_MATRICES = ["BLOSUM62", "BLOSUM45", "BLOSUM50", "BLOSUM80", "BLOSUM90",
                     "PAM30", "PAM70", "PAM250"]  # fmt: skip


@register_function(
    menus="tools/biology/align",
    title="Global Pairwise Alignment",
    command_id="himena-bio:global-pairwise",
)
def global_pairwise_alignment() -> Parametric:
    """Perform a global pairwise alignment of nucleotide sequences."""

    @configure_gui(
        seq0={"types": [Type.SEQS]},
        seq1={"types": [Type.SEQS]},
    )
    def run_global_pairwise_alignment(
        seq0: WidgetDataModel,
        seq1: WidgetDataModel,
        match_score: float = 1.0,
        mismatch_score: float = -0.8,
        gap_score: float = -0.5,
    ) -> WidgetDataModel:
        return _pairwise_impl(
            seq0=seq0,
            seq1=seq1,
            mode="global",
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_score=gap_score,
        )

    return run_global_pairwise_alignment


@register_function(
    menus="tools/biology/align",
    title="Local Pairwise Alignment",
    command_id="himena-bio:local-pairwise",
)
def local_pairwise_alignment() -> Parametric:
    """Perform a local pairwise alignment of nucleotide sequences."""

    @configure_gui(
        seq0={"types": [Type.SEQS]},
        seq1={"types": [Type.SEQS]},
    )
    def run_local_pairwise_alignment(
        seq0: WidgetDataModel,
        seq1: WidgetDataModel,
        match_score: float = 1.0,
        mismatch_score: float = -0.8,
        gap_score: float = -0.5,
    ) -> WidgetDataModel:
        return _pairwise_impl(
            seq0=seq0,
            seq1=seq1,
            mode="local",
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_score=gap_score,
        )

    return run_local_pairwise_alignment


@register_function(
    menus="tools/biology/align",
    title="Protein Pairwise Alignment",
    command_id="himena-bio:protein-pairwise",
)
def protein_pairwise_alignment() -> Parametric:
    """Perform a pairwise alignment of protein sequences using a substitution matrix.

    Default parameters are the same as EMBOSS needle/water (BLOSUM62, gap open -10,
    gap extend -0.5).
    """

    @configure_gui(
        seq0={"types": [Type.PROTEIN]},
        seq1={"types": [Type.PROTEIN]},
        substitution_matrix={"choices": _PROTEIN_MATRICES},
    )
    def run_protein_pairwise_alignment(
        seq0: WidgetDataModel,
        seq1: WidgetDataModel,
        mode: Literal["global", "local"] = "global",
        substitution_matrix: str = "BLOSUM62",
        open_gap_score: float = -10.0,
        extend_gap_score: float = -0.5,
    ) -> WidgetDataModel:
        from Bio.Align import substitution_matrices

        return _pairwise_impl(
            seq0=seq0,
            seq1=seq1,
            mode=mode,
            substitution_matrix=substitution_matrices.load(substitution_matrix),
            open_gap_score=open_gap_score,
            extend_gap_score=extend_gap_score,
        )

    return run_protein_pairwise_alignment


def _pairwise_impl(
    seq0: WidgetDataModel,
    seq1: WidgetDataModel,
    mode: str = "global",
    **kwargs,
) -> WidgetDataModel:
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode=mode, **kwargs)
    seq0_record = current_record(seq0).upper()
    seq1_record = current_record(seq1).upper()
    if len(seq0_record) == 0 or len(seq1_record) == 0:
        raise ValueError("Cannot align empty sequences.")

    alignments = aligner.align(seq0_record, seq1_record)
    return WidgetDataModel(
        value=alignments,
        type=Type.ALIGNMENT,
        title=f"Alignment of {seq0.title} and {seq1.title}",
        extension_default=".aln",
    )
