from himena.plugins import when_command_executed
from himena.workflow import CommandExecution
from himena_bio.consts import Type


def _get_last_id(step: CommandExecution):
    return step.id


(
    when_command_executed(Type.DNA, "himena-bio:pcr").add_command_suggestion(
        "himena-bio:gibson-assembly", defaults={"vec": _get_last_id}
    )
)
(
    when_command_executed(
        Type.DNA, "himena-bio:sanger-sequencing"
    ).add_command_suggestion("himena-bio:local-pairwise")
)
