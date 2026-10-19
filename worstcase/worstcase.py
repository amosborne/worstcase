"""Scalar expression graphs and explicit uncertainty evaluation."""

import ast
import operator
import sys
import warnings
from dataclasses import dataclass
from functools import wraps
from inspect import signature
from itertools import product
from numbers import Integral
from types import MappingProxyType

import numpy as np
from executing import Source

from .unit_adapter import (
    ScalarAdapter,
    dimensionless_value,
    normalize,
    snapshot,
    unit_adapter,
    validate,
)


def _tag(explicit):
    if explicit is not None:
        if not isinstance(explicit, str):
            raise TypeError("tag must be a string or None")
        return explicit
    frame = None
    try:
        frame = sys._getframe(1)
        for _ in range(32):
            if frame is None:
                return ""
            if frame.f_globals.get("__name__") != __name__:
                break
            frame = frame.f_back
        else:
            return ""
        node = Source.executing(frame).node
        parent = getattr(node, "parent", None)
        if isinstance(parent, ast.Assign) and parent.value is node:
            target = parent.targets[0] if len(parent.targets) == 1 else None
        elif (
            isinstance(parent, (ast.AnnAssign, ast.NamedExpr)) and parent.value is node
        ):
            target = parent.target
        else:
            return ""
        return target.id if isinstance(target, ast.Name) else ""
    except Exception:
        # Missing source must never affect a calculation or keep a frame alive.
        return ""
    finally:
        del frame


class _Display:
    def __repr__(self):
        return format(self, "")

    def __format__(self, spec):
        latex = spec.endswith("L")
        spec = spec[:-1] if latex else spec
        spec = spec or ".4g"
        nom, lb, ub = self.nom, self.lb, self.ub
        adapter = unit_adapter(type(nom))

        def render(value):
            return adapter.format(value, spec, latex)

        if latex:
            upper, lower = ub - nom, nom - lb
            if bool(upper == lower):
                return "$" + render(nom) + r" \pm " + render(upper) + "$"
            return (
                "$"
                + render(nom)
                + r"\,{}^{+"
                + render(upper)
                + "}_{-"
                + render(lower)
                + "}$"
            )
        prefix = self.tag + ": " if self.tag else ""
        return prefix + ", ".join(
            f"{render(value)} ({label})"
            for value, label in [(nom, "nom"), (lb, "lb"), (ub, "ub")]
        )

    def _repr_latex_(self):
        return format(self, "L")


@dataclass(frozen=True, eq=False, repr=False)
class Expression:
    """An immutable scalar computation; evaluate bounds explicitly with derive()."""

    _func: object
    _args: tuple
    tag: str = ""
    __array_priority__ = 10000

    def __repr__(self):
        label = self.tag or getattr(self._func, "__name__", "expression")
        return f"<Expression {label}; call .derive() for bounds>"

    def __bool__(self):
        raise TypeError(
            "An uncertain expression has no truth value; use a wrapped function"
        )

    def __array__(self, dtype=None, copy=None):
        raise TypeError("Array-valued uncertainty is unsupported; use @derive")

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__" or kwargs or ufunc not in _UFUNCS:
            return NotImplemented
        return _expression(ufunc, inputs)

    def _binary(self, func, other):
        return _expression(func, (self, other))

    def _reflected(self, func, other):
        return _expression(func, (other, self))

    def __add__(self, other):
        return self._binary(operator.add, other)

    def __radd__(self, other):
        return self._reflected(operator.add, other)

    def __sub__(self, other):
        return self._binary(operator.sub, other)

    def __rsub__(self, other):
        return self._reflected(operator.sub, other)

    def __mul__(self, other):
        return self._binary(operator.mul, other)

    def __rmul__(self, other):
        return self._reflected(operator.mul, other)

    def __truediv__(self, other):
        return self._binary(operator.truediv, other)

    def __rtruediv__(self, other):
        return self._reflected(operator.truediv, other)

    def __pow__(self, other):
        return self._binary(operator.pow, other)

    def __rpow__(self, other):
        return self._reflected(operator.pow, other)

    def __neg__(self):
        return _expression(operator.neg, (self,))

    def __pos__(self):
        return self

    def __abs__(self):
        return _expression(operator.abs, (self,))

    @property
    def nom(self):
        order = _order(self)
        return _evaluate(order, _nominal_inputs(order))

    @property
    def units(self):
        value = self.nom
        return unit_adapter(type(value)).units(value)

    @property
    def u(self):
        return self.units

    def derive(self, method="ev", *, n=1000, seed=None, tag=None):
        """Evaluate a stable result using ev, mc, or rss (see method assumptions)."""
        return _derive(self, method, n, seed, self.tag if tag is None else _tag(tag))

    def ss(self, parameters, method="ev", *, n=1000, seed=None, tag=None):
        """Vary selected primitive parameters/expressions; freeze other inputs."""
        return _derive(
            self, method, n, seed, self.tag if tag is None else _tag(tag), parameters
        )


@dataclass(frozen=True, eq=False, repr=False, init=False)
class Parameter(_Display, Expression):
    """An independent primitive parameter. Equal specifications retain identity."""

    _nom: object
    _lb: object
    _ub: object

    def __init__(self, nom, lb, ub, tag=None):
        validate(nom)
        lb, ub = normalize(lb, nom), normalize(ub, nom)
        if not bool(lb <= nom <= ub):
            raise ValueError("Parameter bounds are out of order")
        object.__setattr__(self, "_func", None)
        object.__setattr__(self, "_args", ())
        object.__setattr__(self, "tag", _tag(tag))
        for name, value in [("_nom", nom), ("_lb", lb), ("_ub", ub)]:
            object.__setattr__(self, name, snapshot(value))

    @property
    def nom(self):
        return snapshot(self._nom)

    @property
    def lb(self):
        return snapshot(self._lb)

    @property
    def ub(self):
        return snapshot(self._ub)

    @staticmethod
    def byrange(nom, lb, ub, tag=None):
        return Parameter(nom, lb, ub, tag)

    @staticmethod
    def constant(value, tag=None):
        return Parameter(value, value, value, tag)

    @staticmethod
    def bytol_rel(nom, rel_tol, tag=None):
        tolerance = dimensionless_value(rel_tol)
        if tolerance < 0:
            raise ValueError("Tolerance must be nonnegative")
        validate(nom)
        delta = abs(nom * tolerance)
        return Parameter(nom, nom - delta, nom + delta, tag)

    @staticmethod
    def bytol_abs(nom, abs_tol, tag=None):
        tolerance = normalize(abs_tol, nom)
        if bool(tolerance < nom * 0):
            raise ValueError("Tolerance must be nonnegative")
        return Parameter(nom, nom - tolerance, nom + tolerance, tag)


@dataclass(frozen=True)
class Derivation:
    """Input assignments; RSS has no joint assignment for either bound."""

    _nom: tuple
    _lb: tuple | None
    _ub: tuple | None

    @staticmethod
    def _mapping(values):
        if values is None:
            return None
        return MappingProxyType({p: snapshot(value) for p, value in values})

    @property
    def nom(self):
        return self._mapping(self._nom)

    @property
    def lb(self):
        return self._mapping(self._lb)

    @property
    def ub(self):
        return self._mapping(self._ub)


@dataclass(frozen=True, eq=False, repr=False)
class Result(_Display):
    """An evaluated snapshot. Reuse the expression to preserve dependencies."""

    _nom: object
    _lb: object
    _ub: object
    tag: str
    method: str
    n_evaluations: int
    seed: int | None
    derivation: Derivation

    @property
    def nom(self):
        return snapshot(self._nom)

    @property
    def lb(self):
        return snapshot(self._lb)

    @property
    def ub(self):
        return snapshot(self._ub)

    @property
    def units(self):
        return unit_adapter(type(self._nom)).units(self._nom)

    @property
    def u(self):
        return self.units


_UFUNCS = {
    np.add,
    np.subtract,
    np.multiply,
    np.divide,
    np.power,
    np.negative,
    np.positive,
    np.absolute,
    np.sqrt,
    np.square,
    np.exp,
    np.log,
    np.log10,
    np.sin,
    np.cos,
    np.tan,
}


_PURE_OPERATIONS = _UFUNCS | {
    operator.add,
    operator.sub,
    operator.mul,
    operator.truediv,
    operator.pow,
    operator.neg,
    operator.abs,
}
_NATIVE_MODULES = {
    "builtins",
    "numpy",
    "pint",
    "astropy",
    "forallpeople",
    "unyt",
    "quantities",
}


def _expression(func, args):
    def node(value):
        return (
            value
            if isinstance(value, Expression)
            else Parameter.constant(value, tag="")
        )

    return Expression(
        func,
        tuple(node(a) for a in args),
        tag=_tag(None),
    )


@dataclass(frozen=True)
class _Slot:
    index: int


def _capture(value, inputs):
    if isinstance(value, Expression):
        inputs.append(value)
        return _Slot(len(inputs) - 1)
    if isinstance(value, (tuple, list)):
        return type(value)(_capture(item, inputs) for item in value)
    if isinstance(value, dict):
        return {k: _capture(v, inputs) for k, v in value.items()}
    return snapshot(value)


def _restore(value, inputs):
    if isinstance(value, _Slot):
        return snapshot(inputs[value.index])
    if isinstance(value, (tuple, list)):
        return type(value)(_restore(item, inputs) for item in value)
    if isinstance(value, dict):
        return {k: _restore(v, inputs) for k, v in value.items()}
    return snapshot(value)


def _has_expression(value):
    if isinstance(value, Expression):
        return True
    if isinstance(value, (tuple, list)):
        return any(_has_expression(item) for item in value)
    if isinstance(value, dict):
        return any(_has_expression(item) for item in value.values())
    return False


def derive(func):
    """Keep concrete calls concrete; build graphs for symbolic scalar calls."""
    sig = signature(func)

    @wraps(func)
    def wrapped(*args, **kwargs):
        binding = sig.bind(*args, **kwargs)
        binding.apply_defaults()
        arguments = (binding.args, binding.kwargs)
        if not _has_expression(arguments):
            return func(*args, **kwargs)
        inputs = []
        template = _capture(arguments, inputs)

        @wraps(func)
        def evaluate(*values):
            positional, keywords = _restore(template, values)
            return func(*positional, **keywords)

        return Expression(evaluate, tuple(inputs), tag=_tag(None))

    return wrapped


def _order(root):
    """Iterative postorder: visit shared nodes once and detect cycles early."""
    ordered, active, done = [], set(), set()
    stack = [(root, False)]
    while stack:
        node, leaving = stack.pop()
        if leaving:
            active.remove(node)
            done.add(node)
            ordered.append(node)
        elif node not in done:
            if node in active:
                raise ValueError("Expression dependencies contain a cycle")
            active.add(node)
            stack.append((node, True))
            children = node._args
            stack.extend((child, False) for child in reversed(children))
    return ordered


def _evaluate(order, assignments):
    values = {}
    for node in order:
        if isinstance(node, Parameter):
            values[node] = snapshot(assignments[node])
        else:
            try:
                operands = [values[a] for a in node._args]
                pure = node._func in _PURE_OPERATIONS and all(
                    type(value).__module__.split(".")[0] in _NATIVE_MODULES
                    for value in operands
                )
                values[node] = validate(
                    node._func(
                        *(operands if pure else map(snapshot, operands)),
                    )
                )
            except Exception as exc:
                label = node.tag or getattr(node._func, "__name__", "expression")
                exc.add_note(f"While evaluating {label!r}")
                raise
    return snapshot(values[order[-1]])


def _nominal_inputs(order):
    assignments = {p: p.nom for p in order if isinstance(p, Parameter)}
    backends, registries = set(), set()
    for value in assignments.values():
        adapter = unit_adapter(type(value))
        if type(adapter) is not ScalarAdapter:
            backends.add(type(adapter))
        registry = getattr(value, "_REGISTRY", None)
        if registry is not None:
            registries.add(id(registry))
    if len(backends) > 1 or len(registries) > 1:
        raise ValueError("An expression must use one unit library and registry")
    return assignments


def _derive(root, method, n, seed, tag, sensitivity=None):
    if method not in ("ev", "mc", "rss"):
        raise ValueError("method must be 'ev', 'mc', or 'rss'")
    if method == "mc":
        if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
            raise ValueError("n must be a positive integer")
        if seed is not None and (
            isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0
        ):
            raise ValueError("seed must be a nonnegative integer or None")
    elif seed is not None or n != 1000:
        raise ValueError("n and seed apply only to Monte Carlo evaluation")
    order = _order(root)
    params = [p for p in order if isinstance(p, Parameter)]
    nominal = _nominal_inputs(order)
    varied = [
        p
        for p in params
        if bool(
            unit_adapter(type(p.nom)).magnitude(p.lb)
            != unit_adapter(type(p.nom)).magnitude(p.ub)
        )
    ]
    if sensitivity is not None:
        selected = (
            [sensitivity] if isinstance(sensitivity, Expression) else list(sensitivity)
        )
        allowed = set()
        for item in selected:
            if not isinstance(item, Expression) or item not in order:
                raise ValueError("Sensitivity inputs must belong to this expression")
            allowed.update(p for p in _order(item) if isinstance(p, Parameter))
        varied = [p for p in varied if p in allowed]

    nom = _evaluate(order, nominal)
    count = 1

    def evaluate(assignments):
        nonlocal count
        try:
            value = normalize(_evaluate(order, assignments), nom)
        except Exception as exc:
            context = ", ".join(
                f"{p.tag or 'input'}={value!s}" for p, value in assignments.items()
            )
            exc.add_note(f"{method.upper()} input assignment: {context}")
            raise
        count += 1
        return value

    lb = ub = nom
    lower_inputs = upper_inputs = nominal.copy()
    if method == "rss":
        lower_squared = upper_squared = None
        asymmetric = False
        for p in varied:
            lo = evaluate({**nominal, p: p.lb}) - nom
            hi = evaluate({**nominal, p: p.ub}) - nom
            zero = nom * 0
            lower = max(-lo, -hi, zero)
            upper = max(lo, hi, zero)
            if bool(max(lower, upper) > (lower + upper) / 2 * 1.05):
                asymmetric = True
            lower_squared = (
                lower**2 if lower_squared is None else lower_squared + lower**2
            )
            upper_squared = (
                upper**2 if upper_squared is None else upper_squared + upper**2
            )
        if asymmetric:
            warnings.warn(
                "RSS is not recommended for asymmetric derived parameters",
                UserWarning,
                stacklevel=3,
            )
        if varied:
            lb, ub = nom - lower_squared**0.5, nom + upper_squared**0.5
        lower_inputs = upper_inputs = None
    else:
        if not varied:
            assignments = ()
        elif method == "ev":
            assignments = (
                {**nominal, **dict(zip(varied, combo))}
                for combo in product(*[(p.lb, p.ub) for p in varied])
            )
        else:
            rng = np.random.default_rng(seed)
            # One uniformly jittered sample per stratum, independently permuted
            # for every primitive: Latin hypercube sampling on uniform ranges.
            columns = [(rng.permutation(n) + rng.random(n)) / n for p in varied]
            ranges = []
            for p in varied:
                lower, upper = p.lb, p.ub
                ranges.append((p, lower, upper - lower))
            assignments = (
                {
                    **nominal,
                    **{
                        p: lower + float(x) * width
                        for (p, lower, width), x in zip(ranges, row)
                    },
                }
                for row in zip(*columns)
            )
        for assignment in assignments:
            value = evaluate(assignment)
            if bool(value < lb):
                lb, lower_inputs = value, assignment
            if bool(value > ub):
                ub, upper_inputs = value, assignment

    def freeze(mapping):
        return (
            None
            if mapping is None
            else tuple((p, snapshot(v)) for p, v in mapping.items())
        )

    validate(lb)
    validate(ub)
    return Result(
        snapshot(nom),
        snapshot(lb),
        snapshot(ub),
        tag,
        method,
        count,
        int(seed) if seed is not None else None,
        Derivation(freeze(nominal), freeze(lower_inputs), freeze(upper_inputs)),
    )
