import logging  # noqa: D100

logger = logging.getLogger(__name__)


def remove_kvl(n):
    """
    Removes Kirchhoff's voltage law (KVL) constraints.

    Function implemented for Kamran's research, and not added to default configs.
    """
    logger.warning("Removing Kirchhoff's Voltage Law constraints")
    n.model.constraints.remove("Kirchhoff-Voltage-Law")
