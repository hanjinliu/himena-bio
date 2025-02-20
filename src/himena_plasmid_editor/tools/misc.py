from himena import WidgetDataModel, StandardType
from himena.plugins import register_function

@register_function(
    menus="tools/biology",
    title="Show Codon Table",
    command_id="himena-plasmid-editor:show-codon-table",
)
def show_codon_table() -> WidgetDataModel:
    from Bio.Seq import CodonTable
    
    return WidgetDataModel(
        value=str(CodonTable.standard_dna_table),
        type=StandardType.TEXT,
        editable=False,
    )
