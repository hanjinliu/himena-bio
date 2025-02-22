from himena import Parametric, WidgetDataModel, StandardType
from himena.consts import MenuId
from himena.plugins import register_function, configure_gui
from himena_bio.consts import Type
from himena_bio.tools._cast import cast_meta, cast_seq_record


@register_function(
    menus=MenuId.FILE_NEW,
    title="New DNA",
    command_id="himena-bio:new-dna",
)
def new_dna() -> WidgetDataModel:
    """Create a new DNA sequence."""
    return WidgetDataModel(value=[""], type=Type.DNA)


@register_function(
    menus="tools/biology",
    title="Show Codon Table",
    command_id="himena-bio:show-codon-table",
)
def show_codon_table() -> WidgetDataModel:
    """Display the standard codon table."""
    from Bio.Seq import CodonTable

    return WidgetDataModel(
        value=str(CodonTable.standard_dna_table),
        type=StandardType.TEXT,
        editable=False,
    )


@register_function(
    menus="tools/biology",
    title="Duplicate selection",
    command_id="himena-bio:duplicate-selection",
)
def duplicate_selection(model: WidgetDataModel) -> Parametric:
    meta = cast_meta(model.metadata)

    @configure_gui(
        current_index={"bind": lambda *_: meta.current_index},
        selection={"bind": lambda *_: meta.selection},
    )
    def run_duplicate(
        current_index: int, selection: tuple[int, int]
    ) -> WidgetDataModel:
        selection_start, selection_end = selection
        original_sequence = cast_seq_record(model.value[current_index])
        new_sequence = original_sequence[selection_start:selection_end]

        return WidgetDataModel(
            value=[new_sequence], type=model.type
        ).with_title_numbering()

    return run_duplicate


@register_function(
    menus="tools/biology",
    title="Duplicate this entry",
    command_id="himena-bio:duplicate-this-entry",
)
def duplicate_this_entry(model: WidgetDataModel) -> Parametric:
    meta = cast_meta(model.metadata)

    @configure_gui(
        current_index={"bind": lambda *_: meta.current_index},
    )
    def run_duplicate(current_index: int) -> WidgetDataModel:
        seq = cast_seq_record(model.value[current_index])
        return WidgetDataModel(value=[seq], type=model.type).with_title_numbering()

    return run_duplicate


@register_function(
    menus="tools/biology",
    title="Reverse complement",
    command_id="himena-bio:reverse-complement",
)
def reverse_complement(model: WidgetDataModel) -> WidgetDataModel:
    out = [cast_seq_record(rec).reverse_complement() for rec in model.value]
    return WidgetDataModel(value=out, type=model.type, title=f"RC of {model.title}")


# TODO: Restriction Digest, Ligation, fetch sequence, etc.


@register_function(
    menus="tools/biology",
    types=[Type.DNA],
    title="PCR",
    command_id="himena-bio:prc",
)
def in_silico_pcr(model: WidgetDataModel) -> Parametric:
    from himena_bio._func import pcr

    def run_pcr(forward: str, reverse: str) -> WidgetDataModel:
        out = []
        for rec in model.value:
            out.append(pcr(rec, forward, reverse))

        return WidgetDataModel(
            value=out, type=model.type, title=f"PCR of {model.title}"
        )

    return run_pcr


@register_function(
    menus="tools/biology",
    title="In-Fusion",
    command_id="himena-bio:in-fusion",
)
def in_silico_in_fusion() -> Parametric:
    from himena_bio._func import in_fusion

    @configure_gui(
        vec={"types": [Type.DNA]},
        insert={"types": [Type.DNA]},
    )
    def run_in_fusion(vec: WidgetDataModel, insert: WidgetDataModel) -> WidgetDataModel:
        out = in_fusion(vec.value[0], insert.value[0])
        return WidgetDataModel(
            value=[out], type=vec.type, title=f"In-Fusion of {vec.title}"
        )

    return run_in_fusion
