from himena import MainWindow, Parametric, WidgetDataModel, StandardType
from himena.widgets import SubWindow
from himena.consts import MenuId
from himena.plugins import register_function, configure_gui
from himena_bio.consts import Type
from himena_bio.tools._cast import cast_meta, cast_seq_record, current_record


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
    types=[Type.SEQS],
    title="Duplicate selection",
    command_id="himena-bio:duplicate-selection",
    group="nucleotide",
)
def duplicate_selection(model: WidgetDataModel) -> Parametric:
    """Duplicate the selected region of the sequence."""
    meta = cast_meta(model.metadata)

    @configure_gui(
        current_index={"bind": lambda *_: meta.current_index},
        selection={"bind": lambda *_: meta.selection},
    )
    def run_duplicate(
        current_index: int, selection: tuple[int, int]
    ) -> WidgetDataModel:
        from himena_bio._utils import slice_seq_record

        selection_start, selection_end = selection
        if selection_start >= selection_end:
            raise ValueError("No region is selected.")
        original_sequence = cast_seq_record(model.value[current_index])
        new_sequence = slice_seq_record(
            original_sequence, slice(selection_start, selection_end)
        )
        new_sequence.annotations["topology"] = "linear"

        return WidgetDataModel(
            value=[new_sequence], type=model.type
        ).with_title_numbering()

    return run_duplicate


@register_function(
    menus="tools/biology",
    types=[Type.SEQS],
    title="Duplicate this entry",
    command_id="himena-bio:duplicate-this-entry",
    group="nucleotide",
)
def duplicate_this_entry(model: WidgetDataModel) -> Parametric:
    """Duplicate the current entry."""
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
    types=[Type.DNA, Type.RNA],
    title="Reverse Complement",
    command_id="himena-bio:reverse-complement",
    group="nucleotide",
)
def reverse_complement(model: WidgetDataModel) -> WidgetDataModel:
    """Reverse complement the sequence."""
    from himena.types import is_subtype

    is_rna = is_subtype(model.type, Type.RNA)
    out = []
    for rec in model.value:
        rec = cast_seq_record(rec)
        rc = rec.reverse_complement(
            id=True, name=True, description=True, annotations=True, dbxrefs=True
        )
        if is_rna:
            rc.seq = rec.seq.reverse_complement_rna()
        out.append(rc)
    return WidgetDataModel(
        value=out,
        type=model.type,
        title=f"RC of {model.title}",
        metadata=model.metadata,
    )


# TODO: Restriction Digest, fetch sequence, etc.


@register_function(
    menus="tools/biology",
    types=[Type.DNA],
    title="PCR",
    command_id="himena-bio:pcr",
    group="nucleotide",
)
def in_silico_pcr(win: SubWindow) -> Parametric:
    """Simulate PCR."""
    from himena_bio._func import pcr

    def run_pcr(forward: str, reverse: str, min_match: int = 15) -> WidgetDataModel:
        # NOTE: output model may change if user ran PCR, and found that the template is
        # not circular, and ran again.
        model = win.to_model()
        out = pcr(current_record(model), forward, reverse, min_match=min_match)
        return WidgetDataModel(
            value=[out], type=model.type, title=f"PCR of {model.title}"
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
            out = gibson_assembly_single(current_record(vec))
        else:
            out = gibson_assembly(current_record(vec), current_record(insert))
        return WidgetDataModel(
            value=[out], type=vec.type, title=f"Gibson of {vec.title}"
        )

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
    title="Self Ligation",
    command_id="himena-bio:self-ligation",
    group="nucleotide",
)
def in_silico_self_ligation(model: WidgetDataModel) -> WidgetDataModel:
    """Simulate self-ligation (circularization) of a linear DNA.

    Both ends are assumed to be ligatable (such as phosphorylated inverse PCR
    products).
    """
    from himena_bio._func import self_ligation

    out = self_ligation(current_record(model))
    return WidgetDataModel(
        value=[out], type=model.type, title=f"Self ligation of {model.title}"
    )


@register_function(
    menus="tools/biology",
    types=[Type.DNA],
    title="Sanger Sequencing",
    command_id="himena-bio:sanger-sequencing",
    group="nucleotide",
)
def in_silico_sequencing(win: SubWindow) -> Parametric:
    """Simulate Sanger sequencing."""
    from himena_bio._func import sequencing

    def run_sequencing(seq: str, length: int = 1000) -> WidgetDataModel:
        model = win.to_model()
        out = sequencing(current_record(model), seq, length=length)
        return WidgetDataModel(
            value=[out], type=model.type, title=f"Sequencing of {model.title}"
        )

    return run_sequencing


@register_function(
    menus="tools/biology",
    types=[Type.PROTEIN],
    title="Protein Properties",
    command_id="himena-bio:protein-properties",
    group="protein",
)
def protein_properties(model: WidgetDataModel) -> WidgetDataModel:
    """Calculate properties of a protein sequence."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    out = []
    for rec in model.value:
        rec = cast_seq_record(rec)
        # stop codons and gaps are not amino acids
        analysis = ProteinAnalysis(
            str(rec.seq).upper().replace("*", "").replace("-", "")
        )
        eps_reduced, eps_cyscys = analysis.molar_extinction_coefficient()
        properties = [
            f"Molecular Weight: {analysis.molecular_weight():.2f}",
            f"Extinction Coefficient (reduced): {eps_reduced}",
            f"Extinction Coefficient (with S-S): {eps_cyscys}",
            f"Isoelectric Point (pI): {analysis.isoelectric_point():.3f}",
            f"Aromaticity: {analysis.aromaticity():.3g}",
            f"Instability Index: {analysis.instability_index():.3g}",
            f"Gravy: {analysis.gravy():.3g}",
        ]
        out.append(f"{rec.name}:\n" + "\n".join(properties))

    return WidgetDataModel(
        value="\n\n".join(out),
        type=StandardType.TEXT,
        title=f"ProtParam of {model.title}",
        editable=False,
    )
