import pytest
from astropy import units as au
from pint import UnitRegistry

from worstcase import derive, param


def test_rss_normalizes_function_output_units():
    u = UnitRegistry()
    x = param.byrange(4 * u.V, 3 * u.V, 5 * u.V)

    def variable_scale(value):
        return value.to(u.mV) if value < 4 * u.V else value.to(u.V)

    result = derive.byrss(x)(variable_scale)()
    assert result.lb.to(u.V).magnitude == pytest.approx(3)
    assert result.ub.to(u.V).magnitude == pytest.approx(5)


def test_astropy_percent_tolerance():
    result = param.bytol_rel(10 * au.V, 10 * au.percent)
    assert result.lb.to_value(au.V) == pytest.approx(9)
    assert result.ub.to_value(au.V) == pytest.approx(11)


def test_negative_nominal_relative_tolerance():
    result = param.bytol_rel(-10, 0.1)
    assert (result.lb, result.ub) == (-11, -9)
