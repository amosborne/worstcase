import forallpeople as si
import pytest
import quantities as pq
import unyt
from astropy import units as u
from pint import UnitRegistry

from worstcase.unit_adapter import (
    AstropyAdapter,
    ForAllPeopleAdapter,
    PintAdapter,
    QuantitiesAdapter,
    ScalarAdapter,
    UnytAdapter,
    unit_adapter,
)

unit = UnitRegistry()
si.environment("default", top_level=False)


@pytest.mark.parametrize(
    "value, ret_cls",
    [
        pytest.param(1, ScalarAdapter),
        pytest.param(1.0, ScalarAdapter),
        pytest.param(1e100, ScalarAdapter),
        pytest.param(1e-100, ScalarAdapter),
        pytest.param(None, ScalarAdapter),
        pytest.param("1", ScalarAdapter),
        pytest.param(1 * unit.A, PintAdapter),
        pytest.param(1 * si.A, ForAllPeopleAdapter),
        pytest.param(1 * u.A, AstropyAdapter),
        pytest.param(1 * unyt.A, UnytAdapter),
        pytest.param(1 * pq.A, QuantitiesAdapter),
    ],
)
def test_unit_adapter_mapping(value, ret_cls):
    assert isinstance(unit_adapter(type(value)), ret_cls)


def test_scalar_adapter():
    val = 2
    adapter = unit_adapter(type(val))
    assert adapter.magnitude(val) == val
    assert adapter.units(val) == 1
    assert adapter.compatible(val, val) is True
    assert adapter.compatible(val, 3) is True
    assert adapter.compatible(val, 3.0) is True
    assert adapter.compatible(val, None) is False


@pytest.mark.parametrize(
    "mag, unit",
    [
        pytest.param(2, 1, id="int"),
        pytest.param(2.0, 1, id="float"),
        pytest.param(2, unit.A, id="pint"),
        pytest.param(2, unit.mA, id="pint-scaled"),
        pytest.param(2, si.A, id="forallpeople"),
        pytest.param(2, u.A, id="astropy"),
        pytest.param(2, u.mA, id="astropy-scaled"),
        pytest.param(2, unyt.A, id="unyt"),
        pytest.param(2, unyt.mA, id="unyt-scaled"),
        pytest.param(2, pq.A, id="python-quantities"),
        pytest.param(2, pq.mA, id="python-quantities-scaled"),
    ],
)
def test_adapter_magnitude_value(mag, unit):
    val = mag * unit
    adapter = unit_adapter(type(val))
    assert adapter.magnitude(val) == mag
    assert adapter.units(val) == unit


@pytest.mark.parametrize(
    "orig_val, targ_mag, targ_unit",
    [
        pytest.param(2 * unit.mA, 0.002, unit.A, id="pint"),
        pytest.param(2 * u.mA, 0.002, u.A, id="astropy"),
        pytest.param(2 * unyt.mA, 0.002, unyt.A, id="unyt"),
        pytest.param(2 * pq.mA, 0.002, pq.A, id="python-quantities"),
    ],
)
def test_adapter_unit_conversion(orig_val, targ_mag, targ_unit):
    adapter = unit_adapter(type(orig_val))
    assert adapter.magnitude(orig_val) != targ_mag
    assert adapter.units(orig_val) != targ_unit

    targ = adapter.convert(orig_val, targ_mag * targ_unit)
    assert adapter.magnitude(targ) == targ_mag
    assert adapter.units(targ) == targ_unit


@pytest.mark.parametrize(
    "val1, val2",
    [
        pytest.param(2, 2.0, id="int-to-float"),
        pytest.param(2 * unit.mA, 2 * unit.A, id="pint-scaled"),
        pytest.param(2 * u.mA, 2 * u.A, id="astropy-scaled"),
        pytest.param(2 * unyt.mA, 2 * unyt.A, id="unyt-scaled"),
        pytest.param(2 * pq.mA, 2 * pq.A, id="python-quantities-scaled"),
    ],
)
def test_adapter_compatible(val1, val2):
    assert unit_adapter(type(val1)).compatible(val1, val2) is True
    assert unit_adapter(type(val2)).compatible(val2, val1) is True


@pytest.mark.parametrize(
    "val1, val2",
    [
        pytest.param(2 * si.A, 2 * si.V, id="forallpeople-diff"),
        pytest.param(2 * si.A, 2, id="forallpeople-scalar"),
        pytest.param(2 * unit.A, 2 * unit.V, id="pint-diff"),
        pytest.param(2 * unit.A, 2, id="pint-scalar"),
        pytest.param(2 * u.A, 2 * u.V, id="astropy-diff"),
        pytest.param(2 * u.A, 2, id="astropy-scalar"),
        pytest.param(2 * unyt.A, 2 * unyt.V, id="unyt-diff"),
        pytest.param(2 * unyt.A, 2, id="unyt-scalar"),
        pytest.param(2 * pq.A, 2 * pq.V, id="python-quantities-diff"),
        pytest.param(2 * pq.A, 2, id="python-quantities-scalar"),
    ],
)
def test_adapter_not_compatible(val1, val2):
    assert unit_adapter(type(val1)).compatible(val1, val2) is False
    assert unit_adapter(type(val2)).compatible(val2, val1) is False


@pytest.mark.parametrize(
    "val",
    [
        pytest.param(4, id="scalar-int"),
        pytest.param(4.0, id="scalar-float"),
        pytest.param(4 * si.A / si.A, id="forallpeople-division"),
        pytest.param(4 * unit.dimensionless, id="pint"),
        pytest.param(4 * unit.A / unit.A, id="pint-division"),
        pytest.param(4 * u.dimensionless_unscaled, id="astropy"),
        pytest.param(4 * u.A / u.A, id="astropy-division"),
        pytest.param(4 * unyt.dimensionless, id="unyt"),
        pytest.param(4 * unyt.A / unyt.A, id="unyt-division"),
        pytest.param(4 * pq.A / pq.A, id="python-quantities-division"),
    ],
)
def test_adapter_dimensionless(val):
    assert unit_adapter(type(val)).is_dimensionless(val) is True


@pytest.mark.parametrize(
    "val",
    [
        pytest.param(4 * si.A, id="forallpeople"),
        pytest.param(4 * unit.A, id="pint"),
        pytest.param(4 * u.A, id="astropy"),
        pytest.param(4 * unyt.A, id="unyt"),
        pytest.param(4 * pq.A, id="python-quantities"),
    ],
)
def test_adapter_not_dimensionless(val):
    assert unit_adapter(type(val)).is_dimensionless(val) is False
