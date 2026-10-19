import inspect
import math
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from worstcase import Expression, derive, param


def test_identity_and_shared_dependencies():
    x = param.bytol_rel(10, 0.1)
    y = param.bytol_rel(10, 0.1)
    assert x is not y and x != y
    assert len({x, y}) == 2
    canceled = (x - x).derive()
    assert (canceled.nom, canceled.lb, canceled.ub) == (0, 0, 0)
    independent = (x - y).derive()
    assert (independent.lb, independent.ub) == (-2, 2)
    a, b = x * 2, x * -2
    result = (a + b).derive()
    assert (result.lb, result.ub) == (0, 0)


def test_evaluation_is_explicit_and_display_is_cheap():
    calls = []

    @derive
    def observed(x):
        calls.append(x)
        return x

    expression = observed(param.byrange(1, 0, 2))
    assert "derive" in repr(expression)
    assert calls == []
    assert not hasattr(expression, "lb")
    assert expression.nom == 1
    assert len(calls) == 1
    result = expression.derive(method="mc", n=10, seed=42)
    count = len(calls)
    assert result.lb == result.lb and result.ub == result.ub
    repr(result)
    format(result, ".3gL")
    assert len(calls) == count
    assert result.n_evaluations == 11


def test_arithmetic_and_reflected_operations():
    x = param.byrange(2, 1, 3)
    assert (3 + x).nom == 5
    assert (3 - x).nom == 1
    assert (3 * x).nom == 6
    assert (3 / x).nom == 1.5
    assert (3**x).nom == 9
    assert (x**2).nom == 4
    assert (-x).nom == -2
    assert (+x).nom == 2
    assert abs(-x).nom == 2
    assert (np.float64(3) * x).nom == 6


@pytest.mark.parametrize("method", ["ev", "mc", "rss"])
def test_constant_expressions(method):
    result = (param.constant(3) + 2).derive(method=method)
    assert (result.nom, result.lb, result.ub) == (5, 5, 5)
    assert result.n_evaluations == 1


def test_ev_includes_nominal_but_does_not_claim_interior_extrema():
    x = param.byrange(0, -1, 1)
    square = (x**2).derive()
    assert (square.nom, square.lb, square.ub) == (0, 0, 1)
    shifted = ((param.byrange(1, 0, 2) - 0.5) ** 2).derive()
    assert shifted.lb == 0.25  # Interior minimum 0 is outside EV's guarantee.


def test_seeded_mc_reproducibility_and_global_rng_isolation():
    expression = param.byrange(0, -1, 2) + param.byrange(2, 1, 3)
    np.random.seed(12)
    expected = np.random.random()
    np.random.seed(12)
    first = expression.derive(method="mc", n=40, seed=42)
    assert np.random.random() == expected
    second = expression.derive(method="mc", n=40, seed=42)
    assert (first.lb, first.ub) == (second.lb, second.ub)
    assert first.derivation.lb == second.derivation.lb
    assert first.seed == 42


def test_latin_hypercube_has_one_sample_per_stratum():
    samples = []

    @derive
    def record(x, y):
        samples.append((x, y))
        return x + y

    x, y = param.byrange(0.5, 0, 1), param.byrange(0.5, 0, 1)
    record(x, y).derive(method="mc", n=20, seed=1)
    for column in zip(*samples[1:]):
        assert sorted(int(x * 20) for x in column) == list(range(20))


def test_mc_one_sample_and_nominal_at_endpoint():
    x = param.byrange(0, 0, 1)
    result = x.derive(method="mc", n=1, seed=42)
    assert result.lb == result.nom == 0
    assert 0 < result.ub < 1
    assert result.derivation.lb[x] == 0


@pytest.mark.parametrize("n", [0, -1, True, 1.5, "10"])
def test_invalid_mc_sample_count(n):
    with pytest.raises(ValueError, match="positive integer"):
        param.constant(1).derive(method="mc", n=n)


@pytest.mark.parametrize("seed", [-1, True, 1.5, "seed"])
def test_invalid_seed(seed):
    with pytest.raises(ValueError, match="seed"):
        param.constant(1).derive(method="mc", seed=seed)


def test_invalid_solver_settings():
    x = param.constant(1)
    with pytest.raises(ValueError, match="method"):
        x.derive(method="unknown")
    with pytest.raises(ValueError, match="only"):
        x.derive(n=1)
    with pytest.raises(ValueError, match="only"):
        x.derive(seed=1)


def test_sensitivity_keeps_dependencies_and_validates_membership():
    x, y = param.byrange(1, 0, 2), param.byrange(2, 0, 4)
    a = x * 2
    expression = a + y
    assert (expression.ss(x).lb, expression.ss(a).ub) == (2, 6)
    assert (expression.ss([]).lb, expression.ss([]).ub) == (4, 4)
    assert expression.ss([x, y]).lb == expression.derive().lb
    with pytest.raises(ValueError, match="belong"):
        expression.ss(param.constant(1))
    with pytest.raises(ValueError, match="belong"):
        expression.ss([1])


def test_signatures_defaults_positional_only_and_errors():
    @derive
    def f(x, /, y=2, *, scale=3):
        return (x + y) * scale

    assert str(inspect.signature(f)) == "(x, /, y=2, *, scale=3)"
    assert f(1) == 9
    assert f(param.byrange(1, 0, 2)).derive().ub == 12
    with pytest.raises(TypeError):
        f(x=1)
    with pytest.raises(TypeError):
        f()

    @derive
    def broken(x):
        raise TypeError("failure inside user function")

    with pytest.raises(TypeError, match="inside user function") as info:
        broken(param.constant(1)).derive()
    assert any("evaluating" in note for note in info.value.__notes__)


def test_cycle_detection_and_deep_graph():
    x = param.constant(1)
    expression = x + 1
    object.__setattr__(expression, "_args", (expression,))
    with pytest.raises(ValueError, match="cycle"):
        expression.derive()
    for _ in range(1200):
        x = x + 1
    assert x.derive().nom == 1201


@pytest.mark.parametrize(
    "func", [np.sqrt, np.exp, np.log, np.log10, np.sin, np.cos, np.tan]
)
def test_numpy_scalar_functions(func):
    x = param.byrange(1, 0.5, 1.5)
    assert func(x).nom == pytest.approx(func(1))


def test_unsupported_implicit_operations_fail_clearly():
    x = param.byrange(1, 0, 2)
    with pytest.raises(TypeError):
        bool(x)
    with pytest.raises(TypeError):
        np.asarray(x)
    with pytest.raises(TypeError):
        np.floor(x)
    with pytest.raises(TypeError):
        np.add(x, 1, where=True)
    with pytest.raises(TypeError):
        x + x.derive()


def test_immutable_nodes_and_independent_provenance():
    x = param.byrange(1, 0, 2)
    with pytest.raises(FrozenInstanceError):
        x.tag = "new"
    result = (x + 1).derive()
    with pytest.raises(FrozenInstanceError):
        result.method = "mc"
    with pytest.raises(TypeError):
        result.derivation.lb[x] = 10
    assert param.constant(2).derive().derivation.nom != result.derivation.nom


def test_rss_analytic_and_asymmetric_results():
    x, y = param.bytol_abs(1, 2), param.bytol_abs(2, 5)
    result = (x * y).derive(method="rss")
    assert result.nom == 2
    assert result.lb == pytest.approx(2 - math.sqrt(41))
    assert result.ub == pytest.approx(2 + math.sqrt(41))
    assert result.derivation.lb is None and result.derivation.ub is None
    with pytest.warns(UserWarning, match="asymmetric"):
        param.byrange(0, 0, 5).derive(method="rss")


def test_name_inference_and_explicit_overrides():
    resistor = param.bytol_rel(10, 0.01)
    assert resistor.tag == "resistor"
    total = resistor + 1
    assert total.tag == "total"
    assert param.constant(1, tag="custom").tag == "custom"
    assert param.constant(1, tag="").tag == ""
    assert eval("param.constant(1)").tag == ""
    with pytest.raises(TypeError, match="tag"):
        param.constant(1, tag=1)


def test_numeric_and_latex_formatting():
    p = param.byrange(1.23456, 1, 2, tag="test")
    assert format(p, ".3g") == "test: 1.23 (nom), 1 (lb), 2 (ub)"
    assert format(p, ".2f") == "test: 1.23 (nom), 1.00 (lb), 2.00 (ub)"
    assert format(p, ".3gL") == r"$1.23\,{}^{+0.765}_{-0.235}$"
    assert format(param.bytol_abs(1, 0.25), ".2fL") == r"$1.00 \pm 0.25$"
    assert p._repr_latex_().startswith("$")
    assert p.derive()._repr_latex_().endswith("$")


def test_public_expression_type():
    assert isinstance(param.constant(1) + 1, Expression)


def test_wrapped_static_configuration_and_nested_scalar_inputs():
    @derive
    def configured(values, mode="double", *, options=None):
        return sum(values) * (2 if mode == "double" else 1) + options["offset"]

    x = param.byrange(2, 1, 3)
    expression = configured([x, x], options={"offset": 1})
    result = expression.derive()
    assert (result.nom, result.lb, result.ub) == (9, 5, 13)
    assert result.n_evaluations == 3


def test_wrapped_callback_mutation_cannot_change_expression_inputs():
    @derive
    def consume(value, items):
        return value + items.pop()

    items = [1]
    expression = consume(param.byrange(2, 1, 3), items)
    items.append(100)
    result = expression.derive()
    assert (result.nom, result.lb, result.ub) == (3, 2, 4)
