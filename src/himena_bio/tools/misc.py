from himena import MainWindow, Parametric, WidgetDataModel, StandardType
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
    group="others",
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
    group="nucleotide",
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
    group="nucleotide",
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
    title="Reverse Complement",
    command_id="himena-bio:reverse-complement",
    group="nucleotide",
)
def reverse_complement(model: WidgetDataModel) -> WidgetDataModel:
    """Reverse complement the sequence."""
    out = [cast_seq_record(rec).reverse_complement() for rec in model.value]
    return WidgetDataModel(value=out, type=model.type, title=f"RC of {model.title}")


# TODO: Restriction Digest, Ligation, fetch sequence, etc.


@register_function(
    menus="tools/biology",
    types=[Type.DNA],
    title="PCR",
    command_id="himena-bio:pcr",
    group="nucleotide",
)
def in_silico_pcr(model: WidgetDataModel) -> Parametric:
    """Simulate PCR."""
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
    title="Gibson Assembly",
    command_id="himena-bio:gibson-assembly",
    group="nucleotide",
)
def in_silico_gibson_assembly() -> Parametric:
    """Simulate cloning by Gibson assembly."""
    from himena_bio._func import gibson_assembly, gibson_assembly_single

    @configure_gui(
        vec={"types": [Type.DNA]},
        insert={"types": [Type.DNA]},
    )
    def run_gibson(
        vec: WidgetDataModel,
        insert: WidgetDataModel | None = None,
    ) -> WidgetDataModel:
        if insert is None:
            out = gibson_assembly_single(vec.value[0])
        else:
            out = gibson_assembly(vec.value[0], insert.value[0])
        return WidgetDataModel(value=[out], type=vec.type)

    return run_gibson


@register_function(
    menus=[],
    types=[Type.DNA],
    title="Gibson Assembly Using This ...",
    command_id="himena-bio:gibson-assembly-this",
    group="nucleotide",
)
def in_silico_gibson_assembly_this(model: WidgetDataModel, ui: MainWindow):
    """Simulate cloning by Gibson assembly."""
    return ui.exec_action("himena-bio:gibson-assembly", with_defaults={"vec": model})


@register_function(
    menus="tools/biology",
    types=[Type.DNA],
    title="Sanger Sequencing",
    command_id="himena-bio:sanger-sequencing",
    group="nucleotide",
)
def in_silico_sequencing(model: WidgetDataModel) -> Parametric:
    """Simulate Sanger sequencing."""
    from himena_bio._func import sequencing

    def run_sequencing(seq: str) -> WidgetDataModel:
        out = []
        for rec in model.value:
            out.append(sequencing(rec, seq))

        return WidgetDataModel(
            value=out, type=model.type, title=f"Sequencing of {model.title}"
        )

    return run_sequencing
