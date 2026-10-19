import math

import numpy as np
import pytest

from worstcase import derive, param
from worstcase.unit_adapter import dimensionless_value, normalize, unit_adapter


def relative(value, reference):
    return dimensionless_value(value / reference)


@pytest.mark.parametrize("method", ["ev", "mc", "rss"])
def test_all_backends_and_scaled_bounds(backend, method):
    v, mv = backend.V, backend.mV
    x = param.byrange(4 * v, 3900 * mv, 4.1 * v)
    y = param.bytol_abs(2 * v, 100 * mv)
    result = (x + y).derive(method=method, **({"seed": 42} if method == "mc" else {}))
    assert relative(result.nom, 6 * v) == pytest.approx(1)
    delta = math.sqrt(0.02) if method == "rss" else 0.2
    tolerance = 0.01 if method == "mc" else 1e-9
    assert relative(result.lb, (6 - delta) * v) == pytest.approx(1, abs=tolerance)
    assert relative(result.ub, (6 + delta) * v) == pytest.approx(1, abs=tolerance)


def test_tolerances_and_percentages(backend):
    v = backend.V
    negative = param.bytol_rel(-10 * v, 10 * backend.percent)
    assert relative(negative.lb, -11 * v) == pytest.approx(1)
    assert relative(negative.ub, -9 * v) == pytest.approx(1)
    with pytest.raises(ValueError, match="nonnegative"):
        param.bytol_rel(10 * v, -0.1)
    with pytest.raises(ValueError, match="nonnegative"):
        param.bytol_abs(10 * v, -1 * v)
    with pytest.raises(ValueError, match="order"):
        param.byrange(1 * v, 2 * v, 3 * v)


def test_incompatible_dimensions_and_bare_numbers(backend):
    if backend.name == "plain":
        return
    with pytest.raises(ValueError):
        param.byrange(1 * backend.V, 0.5, 1.5)
    with pytest.raises(ValueError):
        param.bytol_abs(1 * backend.V, 1 * backend.A)
    with pytest.raises(ValueError, match="dimensionless"):
        param.bytol_rel(1 * backend.V, 0.1 * backend.A)


def test_quantity_constants_and_provenance(backend):
    v = backend.V
    x = param.byrange(2 * v, 1 * v, 3 * v)
    result = (param.constant(5 * v) * x / (x + x)).derive()
    assert relative(result.nom, 2.5 * v) == pytest.approx(1)
    assert relative(result.lb, 2.5 * v) == pytest.approx(1)
    result = (x + 2 * v).derive()
    assert relative(result.lb, 3 * v) == pytest.approx(1)
    assert relative(result.derivation.lb[x], 1 * v) == pytest.approx(1)
    assert bool(result.units == result.u)


def test_output_normalization_rss_regression(backend):
    if backend.name in ("plain", "forallpeople"):
        return
    v, mv = backend.V, backend.mV
    x = param.byrange(4 * v, 3 * v, 5 * v)

    @derive
    def scaled(value):
        return normalize(value, 1 * mv if bool(value < 4 * v) else 1 * v)

    result = scaled(x).derive(method="rss")
    assert relative(result.lb, 3 * v) == pytest.approx(1)
    assert relative(result.ub, 5 * v) == pytest.approx(1)


def test_quantity_storage_is_owned(backend):
    if backend.name in ("plain", "forallpeople"):
        return
    value = 2 * backend.V
    p = param.byrange(value, 1 * backend.V, 3 * backend.V)
    result = (p * 2).derive()
    # All these mutable quantity implementations expose an ndarray magnitude.
    adapter = unit_adapter(type(value))
    magnitude = adapter.magnitude(value)
    if isinstance(magnitude, np.ndarray):
        magnitude[...] = 999
    exposed = result.nom
    if isinstance(adapter.magnitude(exposed), np.ndarray):
        adapter.magnitude(exposed)[...] = -999
    assert relative(p.nom, 2 * backend.V) == pytest.approx(1)
    assert relative(result.nom, 4 * backend.V) == pytest.approx(1)


def test_backend_text_and_latex(backend):
    result = param.bytol_rel(10 * backend.V, 0.1, tag="voltage").derive()
    assert "nom" in format(result, ".3g")
    latex = format(result, ".3gL")
    assert latex.startswith("$") and latex.endswith("$")
    assert r"\pm" in latex


@pytest.mark.parametrize(
    "value", [True, np.bool_(True), "1", None, [], 1j, np.array([1])]
)
def test_unsupported_scalar_types(value):
    with pytest.raises(TypeError):
        param.constant(value)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_inputs(value):
    with pytest.raises(ValueError, match="finite"):
        param.constant(value)


def test_subclass_dispatch_and_dimensionless_coercion():
    from astropy import units as u

    class Voltage(u.Quantity):
        pass

    q = Voltage(1, u.V)
    assert type(unit_adapter(type(q))).__name__ == "AstropyAdapter"
    p = param.byrange(5 * u.dimensionless_unscaled, 4, 6)
    assert p.lb == 4 * u.dimensionless_unscaled
    assert dimensionless_value(10 * u.percent) == pytest.approx(0.1)


def test_cross_library_and_registry_rejection():
    from astropy import units as au
    from pint import UnitRegistry

    u, other = UnitRegistry(), UnitRegistry()
    with pytest.raises(ValueError):
        param.byrange(2 * u.V, 1 * au.V, 3 * au.V)
    with pytest.raises(ValueError):
        param.byrange(2 * u.V, 1 * other.V, 3 * other.V)
    expression = param.constant(1 * u.V) + param.constant(1 * other.V)
    with pytest.raises(ValueError, match="one unit library and registry"):
        expression.nom


def test_forallpeople_small_tolerances_are_not_frozen():
    import forallpeople as fp

    fp.environment("default", top_level=False)
    x = param.bytol_rel(1 * fp.V, 1e-12)
    result = x.derive()
    assert result.n_evaluations == 3
    assert result.lb.value < 1 < result.ub.value
