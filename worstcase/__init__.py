from .worstcase import Expression
from .worstcase import Parameter as param
from .worstcase import Result, derive

__all__ = ["param", "derive", "Expression", "Result"]


def __getattr__(name):
    """
    Provide a more informative error message for importing removed methods.
    """
    if name == "unit":
        raise ImportError("worstcase no longer supplies a unit registry.")
    raise AttributeError(f"Module {__name__} has no attribute {name}")
