import functools
import numbers
from operator import attrgetter


class ScalarAdapter:
    """Common set of methods to adapt different unit libraries."""

    magnitude = staticmethod(lambda q: q)
    """Return the magnitude of a value."""

    units = staticmethod(lambda q: 1)
    """Return the units of a value.

    This function returns the element that can be multiplied to a scalar
        to get the dimensioned value. This does not return a string
        representation of the value's units.
    """

    @staticmethod
    def compatible(a, b):
        """Check if two values are unit-compatible"""
        if isinstance(a, numbers.Real) and isinstance(b, numbers.Real):
            return True
        if not isinstance(a, type(b)):
            return False

        try:
            a + b
            return True
        except (ValueError, TypeError):
            return False

    @staticmethod
    def convert(a, b):
        """Converts `a` to the units of `b`."""
        return a

    @staticmethod
    def format_str(sigfig):
        """Return a valid formatting string for the unit type"""
        return "0.{sigfig}g".format(sigfig=sigfig)

    @classmethod
    def is_dimensionless(cls, value):
        """Check if a value is dimensionless"""
        return True


class PintAdapter(ScalarAdapter):
    """Adapter for Pint Quantities"""

    units = staticmethod(attrgetter("units"))
    magnitude = staticmethod(attrgetter("magnitude"))

    @staticmethod
    def compatible(a, b):
        if not isinstance(a, type(b)):
            return False
        return a.is_compatible_with(b)

    @staticmethod
    def convert(a, b):
        return a.to(b.units)

    @staticmethod
    def format_str(sigfig):
        return "0.{sigfig}G#~P".format(sigfig=sigfig)

    @classmethod
    def is_dimensionless(cls, value):
        return value.dimensionless


class AstropyAdapter(ScalarAdapter):
    """Adapter for astropy Quantity instance"""

    units = staticmethod(attrgetter("unit"))
    magnitude = staticmethod(attrgetter("value"))

    @staticmethod
    def compatible(a, b):
        if not isinstance(a, type(b)):
            return False
        return a.unit.is_equivalent(b.unit)

    @classmethod
    def convert(cls, a, b):
        new_unit = cls.units(b)
        return a.to(new_unit)

    @classmethod
    def is_dimensionless(cls, value):
        return str(cls.units(value)) == ""


class ForAllPeopleAdapter(ScalarAdapter):
    """Adapter for forallpeople Physical instance"""

    units = staticmethod(lambda v: v.split()[1])
    magnitude = staticmethod(attrgetter("value"))

    @classmethod
    def is_dimensionless(cls, value):
        return False


class UnytAdapter(ScalarAdapter):
    """Adapter for unyt instance"""

    units = staticmethod(attrgetter("units"))
    magnitude = staticmethod(attrgetter("value"))

    @staticmethod
    def compatible(a, b):
        if not isinstance(a, type(b)):
            return False
        return a.units.same_dimensions_as(b.units)

    @classmethod
    def convert(cls, a, b):
        new_unit = cls.units(b)
        return a.to(new_unit)

    @classmethod
    def is_dimensionless(cls, value):
        return value.units.is_dimensionless


class QuantitiesAdapter(ScalarAdapter):
    """Adapter for python-quantities Quantity instance"""

    units = staticmethod(attrgetter("units"))
    magnitude = staticmethod(attrgetter("magnitude"))

    @classmethod
    def convert(cls, a, b):
        new_unit = cls.units(b)
        return a.rescale(new_unit)

    @classmethod
    def is_dimensionless(cls, value):
        return "dimensionless" in str(cls.units(value))


_ADAPTERS = {
    "pint": PintAdapter(),
    "astropy": AstropyAdapter(),
    "unyt": UnytAdapter(),
    "quantities": QuantitiesAdapter(),
    "forallpeople": ForAllPeopleAdapter(),
}

_DEFAULT_ADAPTER = ScalarAdapter()


@functools.cache
def unit_adapter(value_type):
    """Map a unit library to an Adapter.

    For a value `v`, get the corresponding unit library adapter with
        `unit_adapter(type(v))`.

    Arguments:
        value_type: Hashable object which has a valid `__module__`.

    Returns:
        Adapter for the provided unit library.
        Returns ScalarAdapter if unit library not found.
    """
    unit_library = value_type.__module__.split(".")[0]
    return _ADAPTERS.get(unit_library, _DEFAULT_ADAPTER)
