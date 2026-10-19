"""Optional quantity integration. Numerical methods operate on quantities directly."""

import functools
import math
import numbers
from copy import copy
from operator import attrgetter

import numpy as np


class ScalarAdapter:
    magnitude = staticmethod(lambda value: value)
    units = staticmethod(lambda value: 1)
    convert = staticmethod(lambda value, reference: value)
    is_dimensionless = staticmethod(lambda value: True)
    dimensionless_value = staticmethod(lambda value: value)

    @staticmethod
    def compatible(a, b):
        return isinstance(a, numbers.Real) and isinstance(b, numbers.Real)

    @classmethod
    def format(cls, value, spec, latex=False):
        return format(value, spec)


class PintAdapter(ScalarAdapter):
    magnitude = staticmethod(attrgetter("magnitude"))
    units = staticmethod(attrgetter("units"))
    is_dimensionless = staticmethod(lambda value: bool(value.dimensionless))
    dimensionless_value = staticmethod(lambda value: value.to("").magnitude)

    @staticmethod
    def compatible(a, b):
        return getattr(a, "_REGISTRY", None) is getattr(
            b, "_REGISTRY", None
        ) and a.is_compatible_with(b)

    @staticmethod
    def convert(value, reference):
        return value.to(reference.units)

    @staticmethod
    def format(value, spec, latex=False):
        return format(value, spec + ("~L" if latex else "~P"))


class AstropyAdapter(ScalarAdapter):
    magnitude = staticmethod(attrgetter("value"))
    units = staticmethod(attrgetter("unit"))
    is_dimensionless = staticmethod(lambda value: value.unit.is_equivalent(1))
    dimensionless_value = staticmethod(lambda value: value.to_value(""))
    compatible = staticmethod(lambda a, b: a.unit.is_equivalent(b.unit))
    convert = staticmethod(lambda value, reference: value.to(reference.unit))

    @staticmethod
    def format(value, spec, latex=False):
        text = format(value.value, spec)
        unit = (
            value.unit.to_string("latex_inline").strip("$")
            if latex
            else str(value.unit)
        )
        return text + ((r"\," if latex else " ") + unit if unit else "")


class ForAllPeopleAdapter(ScalarAdapter):
    magnitude = staticmethod(attrgetter("value"))
    units = staticmethod(lambda value: value.split()[1])
    is_dimensionless = staticmethod(lambda value: not any(value.dimensions))
    dimensionless_value = staticmethod(lambda value: value.value / value.factor)
    compatible = staticmethod(lambda a, b: a.dimensions == b.dimensions)

    @staticmethod
    def format(value, spec, latex=False):
        # Physical's formatter handles engineering prefixes and its own units.
        return format(value, spec + ("L" if latex else ""))


class UnytAdapter(ScalarAdapter):
    magnitude = staticmethod(attrgetter("value"))
    units = staticmethod(attrgetter("units"))
    is_dimensionless = staticmethod(lambda value: bool(value.units.is_dimensionless))
    dimensionless_value = staticmethod(lambda value: value.to_value("dimensionless"))
    compatible = staticmethod(lambda a, b: a.units.same_dimensions_as(b.units))
    convert = staticmethod(lambda value, reference: value.to(reference.units))

    @staticmethod
    def format(value, spec, latex=False):
        text = format(value.value.item(), spec)
        unit = value.units.latex_repr if latex else str(value.units)
        return text + (r"\," if latex else " ") + unit


class QuantitiesAdapter(ScalarAdapter):
    magnitude = staticmethod(attrgetter("magnitude"))
    units = staticmethod(attrgetter("units"))
    is_dimensionless = staticmethod(lambda value: not value.simplified.dimensionality)
    dimensionless_value = staticmethod(lambda value: value.simplified.magnitude)
    compatible = staticmethod(
        lambda a, b: a.simplified.dimensionality == b.simplified.dimensionality
    )
    convert = staticmethod(lambda value, reference: value.rescale(reference.units))

    @staticmethod
    def format(value, spec, latex=False):
        text = format(value.magnitude.item(), spec)
        unit = (
            value.dimensionality.latex.strip("$")
            if latex
            else str(value.dimensionality)
        )
        return text + (r"\," if latex else " ") + unit


_ADAPTERS = {
    "pint": PintAdapter(),
    "astropy": AstropyAdapter(),
    "forallpeople": ForAllPeopleAdapter(),
    "unyt": UnytAdapter(),
    "quantities": QuantitiesAdapter(),
}
_SCALAR = ScalarAdapter()


@functools.cache
def unit_adapter(value_type):
    """Resolve supported quantity subclasses without importing optional packages."""
    for cls in value_type.__mro__:
        adapter = _ADAPTERS.get(cls.__module__.split(".")[0])
        if adapter is not None:
            return adapter
    if issubclass(value_type, numbers.Real) and not issubclass(
        value_type, (bool, np.bool_)
    ):
        return _SCALAR
    raise TypeError(f"Unsupported scalar or quantity type: {value_type.__name__}")


def validate(value):
    """Require a finite, real scalar in a supported backend."""
    adapter = unit_adapter(type(value))
    magnitude = np.asarray(adapter.magnitude(value))
    if magnitude.ndim != 0 or magnitude.dtype.kind not in "iuf":
        raise TypeError("Expected a real scalar, not an array or complex value")
    if not math.isfinite(float(magnitude)):
        raise ValueError("Values must be finite")
    return value


def dimensionless_value(value):
    validate(value)
    adapter = unit_adapter(type(value))
    if not adapter.is_dimensionless(value):
        raise ValueError("Relative tolerance must be dimensionless")
    return float(adapter.dimensionless_value(value))


def normalize(value, reference):
    """Validate and convert a scalar to the reference's backend and units."""
    validate(value)
    validate(reference)
    adapter = unit_adapter(type(value))
    target = unit_adapter(type(reference))
    if adapter is not target:
        if adapter is _SCALAR and target.is_dimensionless(reference):
            units = target.units(reference)
            return target.convert(value * (units / units), reference)
        if target is _SCALAR and adapter.is_dimensionless(value):
            return dimensionless_value(value)
        raise ValueError("Values must use compatible units and the same unit library")
    try:
        if not adapter.compatible(value, reference):
            raise ValueError("Values have incompatible units or registries")
        return adapter.convert(value, reference)
    except (TypeError, AttributeError) as exc:
        raise ValueError("Values have incompatible units or registries") from exc


def snapshot(value):
    """Own quantity storage; caller mutation cannot change stored results."""
    return copy(value)
