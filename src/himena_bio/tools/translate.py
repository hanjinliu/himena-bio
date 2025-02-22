from himena import WidgetDataModel, Parametric
from himena.plugins import register_function, configure_gui
from himena_bio.consts import Type
from himena_bio.tools._cast import cast_meta, cast_seq_record
from himena_bio._utils import topology


@register_function(
    menus="tools/biology",
    types=[Type.DNA, Type.RNA],
    title="Translate",
    keybindings=["Alt+T"],
    command_id="himena-bio:translate",
)
def translate(model: WidgetDataModel) -> Parametric:
    """Translate a DNA or RNA sequence to protein."""
    meta = cast_meta(model.metadata)

    @configure_gui(
        index={"bind": meta.current_index},
        selection={"bind": meta.selection},
    )
    def run_translate(index: int, selection: tuple[int, int]) -> WidgetDataModel:
        record = cast_seq_record(model.value[index])
        start, end = selection
        translations = record.seq[start:end].translate()
        return WidgetDataModel(value=[str(translations)], type=Type.PROTEIN)

    return run_translate


@register_function(
    menus="tools/biology",
    types=[Type.DNA, Type.RNA],
    title="Translate until stop codon",
    keybindings=["Alt+Shift+T"],
    command_id="himena-bio:translate-until-stop",
)
def translate_until_stop(model: WidgetDataModel) -> Parametric:
    """Translate a DNA or RNA sequence from the selection until a stop codon is encountered."""
    meta = cast_meta(model.metadata)

    @configure_gui(
        index={"bind": meta.current_index},
        start={"bind": meta.selection[0]},
    )
    def run_translate(index: int, start: int) -> WidgetDataModel:
        record = cast_seq_record(model.value[index])
        _topo = topology(record)
        if _topo == "linear":
            seq_ref = record.seq[start:]
        elif _topo == "circular":
            seq_ref = record.seq[start:] + record.seq[:start]
        else:
            raise ValueError(f"Invalid topology: {_topo!r}")
        translations = seq_ref.translate(to_stop=True)
        return WidgetDataModel(
            value=[str(translations)],
            type=Type.PROTEIN,
        )

    return run_translate
