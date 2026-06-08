from .worstcase import Derivative as derive
from .worstcase import Parameter as param

__all__ = ["param", "derive"]


def __getattr__(name):
    """
    Provide a more informative error message for importing removed methods.
    """
    if name == "unit":
        raise ImportError("worstcase no longer supplies a unit registry.")
    raise AttributeError(f"Module {__name__} has no attribute {name}")
