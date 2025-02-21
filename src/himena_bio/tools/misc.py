from himena import WidgetDataModel, StandardType, Parametric
from himena.consts import MenuId
from himena.plugins import register_function, configure_gui
from himena_bio.consts import SeqMeta, Type


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
    types=[Type.DNA, Type.RNA],
    title="Translate",
    command_id="himena-bio:translate",
)
def translate(model: WidgetDataModel) -> Parametric:
    """Translate a DNA sequence to protein."""
    from Bio.SeqRecord import SeqRecord

    meta = model.metadata
    if not isinstance(meta, SeqMeta):
        raise ValueError("Invalid metadata")

    @configure_gui(
        index={"bind": meta.current_index},
        selection={"bind": meta.selection},
    )
    def run_translate(index: int, selection: tuple[int, int]):
        record = model.value[index]
        assert isinstance(record, SeqRecord)
        start, end = selection
        translations = record.seq[start:end].translate()
        return WidgetDataModel(value=[str(translations)], type=Type.PROTEIN)

    return run_translate


# TODO: PCR, InFusion, alignment, Restriction Digest, Ligation, fetch sequence, etc.
