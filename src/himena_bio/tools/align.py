from himena import Parametric, WidgetDataModel
from himena.widgets import TabArea
from himena.plugins import register_function, configure_gui
from himena_bio.consts import Type
from himena_bio.tools._cast import cast_seq_record


@register_function(
    menus="tools/biology/align",
    title="Global Pairwise Alignment",
    command_id="himena-bio:global-pairwise",
)
def global_pairwise_alignment(tab: TabArea) -> Parametric:
    """Perform a global pairwise alignment."""

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
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_score=gap_score,
            mode="global",
        )

    return run_global_pairwise_alignment


@register_function(
    menus="tools/biology/align",
    title="Local Pairwise Alignment",
    command_id="himena-bio:local-pairwise",
)
def local_pairwise_alignment() -> Parametric:
    """Perform a local pairwise alignment."""

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
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_score=gap_score,
            mode="local",
        )

    return run_local_pairwise_alignment


def _pairwise_impl(
    seq0: WidgetDataModel,
    seq1: WidgetDataModel,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_score: float = -0.5,
    mode: str = "global",
) -> WidgetDataModel:
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score,
        mode=mode,
    )
    if len(seq0.value) != 1 or len(seq1.value) != 1:
        raise ValueError("Both sequences must be single entry.")
    seq0_record = cast_seq_record(seq0.value[0])
    seq1_record = cast_seq_record(seq1.value[0])

    alignments = aligner.align(seq0_record, seq1_record)
    return WidgetDataModel(
        value=alignments,
        type=Type.ALIGNMENT,
        title="alignment.aln",
        extension_default=".aln",
    )
