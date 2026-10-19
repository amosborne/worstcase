# worstcase

Unit-aware scalar expressions for corner analysis (EV), uniform Monte Carlo
sampling (MC), and endpoint root-sum-square analysis (RSS).

Version 1.0 is a breaking redesign, currently a release candidate. See
[migration instructions](MIGRATION.md) and the [current-sensing notebook](example.ipynb).

## Installation

```sh
pip install .
```

Install this unpublished candidate from a checkout. After publication, use
`pip install worstcase` for stable releases.

Requires Python 3.11 or later. Unit libraries are optional: plain real numbers,
Pint, forallpeople, Astropy, unyt, and python-quantities are supported. Install
your preferred library separately. Use one library/registry per calculation.

## Expressions and explicit evaluation

```python
from pint import UnitRegistry
from worstcase import param, derive

u = UnitRegistry()
R1 = param.bytol_rel(10 * u.kohm, 0.01)
R2 = param.bytol_rel(10 * u.kohm, 0.01)
vin = param.constant(5 * u.V)

vout = vin * R2 / (R1 + R2)
result = vout.derive()
print(f"{result:.4g}")
assert abs(result.nom.to("V").magnitude - 2.5) < 1e-12
```

Create independent primitives with `param.byrange(nom, lb, ub)`,
`param.bytol_rel(nom, rel_tol)` (0.01 means 1%), or
`param.bytol_abs(nom, abs_tol)`. Relative tolerance uses the absolute magnitude of
nominal, including for negative values. Tolerances must be nonnegative. Compatible
bounds are converted to nominal's units; percentages are valid relative tolerances.

Repeated use of a parameter preserves its identity: `R1 - R1` is exactly zero.
Two separately constructed resistors with identical specifications vary independently.
Expressions and parameters are immutable, and input quantities are copied at creation.
Unit registry metadata is shared; configure registries before constructing analyses.

Expressions support `+`, `-`, `*`, `/`, `**`, unary signs, `abs`, and scalar NumPy
arithmetic, `sqrt`, `exp`, `log`, `log10`, `sin`, `cos`, and `tan`. Their `.nom`
property computes nominal only. Expression display does not run a solver; bounds
require `.derive()`. Results are snapshots: reading `.nom`, `.lb`, `.ub`, formatting,
or displaying them never reruns an analysis.

```python
mc = vout.derive(method="mc", n=1000, seed=42)
rss = vout.derive(method="rss")
print(mc.n_evaluations, mc.seed)
print(f"{mc:.3gL}")
sensitivity = vout.ss(R1)
```

`.ss(parameter_or_expression)` varies the selected primitive inputs while holding
all other inputs at nominal. Pass a list for a group, or `[]` to freeze everything.
`method`, `n`, and `seed` work as in `.derive()`. `n` and `seed` apply only to MC.

EV and MC results expose `.derivation.nom`, `.derivation.lb`, and `.derivation.ub`:
read-only maps from primitive objects to input assignments. Ties retain the first
assignment encountered (nominal is evaluated first). RSS has no single joint input
assignment for its combined bounds, so its lower/upper derivations are `None`.
`n_evaluations` counts complete root evaluations, including nominal. Fixed inputs
are excluded from corner enumeration and sampling; constant expressions evaluate once.

## Functions and units

```python
@derive
def gain_with_offset(gain, /, offset=0):
    return gain * 2 + offset

x = param.byrange(2, 1, 3)
answer = gain_with_offset(x, offset=1).derive()
assert (answer.nom, answer.lb, answer.ub) == (5, 3, 7)
assert gain_with_offset(2, offset=1) == 5
```

Wrapped functions receive concrete values during evaluation and return a real scalar
or supported scalar quantity. They can use arrays, conditionals, and numerical
solvers internally. Keep callbacks pure and deterministic. Their exceptions propagate
with evaluation context. Signature/default/positional-only behavior is preserved.

An expression wraps scalar constants automatically when it receives an operation:
`R1 * (5 * u.V)` works. For portable quantity-left arithmetic, wrap the quantity
explicitly as `param.constant(5 * u.V)`: some unit libraries raise or absorb an
expression before Python calls its reflected operator. We do not patch those libraries.
For Pint, expression types register through its documented upcast-type extension
when Pint is loaded. This prevents wrapping an uncertain magnitude and accidentally
discarding its outer units; no library methods are replaced.

Array-valued uncertain inputs/outputs, implicit conversion to arrays, symbolic truth
tests, and unsupported NumPy operations fail explicitly. Use a wrapped scalar-output
function for a quantity of interest. Reuse an expression to preserve dependencies;
results cannot be inserted into expression arithmetic. Constructing a new parameter
from result bounds intentionally creates a new independent uncertainty.

## Display and names

Assignment names are inferred on a best-effort basis, including in ordinary notebooks.
`tag="custom label"` overrides inference; `tag=""` suppresses a label. Source/AST
availability, factory calls, and dynamic code can prevent inference; the fallback
is an empty label. Names never affect identity or calculations.

Use standard numeric format specifications such as `.3g`, `.2f`, or `.2e`.
Append `L` for a LaTeX uncertainty representation (`.3gL`); notebook results also
implement `_repr_latex_`. An omitted specification uses `.4g`. Formatting operates
on result values and never changes their precision or units.

## What the methods mean

- **EV** evaluates every corner of the input box plus nominal. With `k` varying
  primitives, this costs `2**k + 1` evaluations. It finds extrema for functions whose
  extrema lie at those corners; it does not find arbitrary interior extrema or
  prove a global bound. For example, `(x - 0.5)**2` with nominal 1 over `[0, 2]`
  reports lower bound 0.25 even though its true minimum is 0.
- **MC** uses independent uniform ranges with Latin hypercube sampling. The result
  is the observed envelope of `n` samples and nominal, not a confidence interval or
  guaranteed bound. A seed makes runs repeatable for the same expression structure
  and dependency versions, without changing NumPy's global random state.
- **RSS** varies one primitive at a time to both endpoints, then combines positive
  and negative deviations in quadrature. It assumes independent contributions and
  is an approximation for nonlinear functions. It warns about asymmetric endpoint
  responses; absence of a warning does not establish linearity. These ranges do not
  specify standard deviations or statistical confidence levels.

All methods operate on the same original primitives, without treating intermediate
expressions as newly independent inputs. Every evaluated value must remain finite
and dimensionally compatible with the nominal output. The libraries govern affine
temperature and logarithmic unit arithmetic; use linear units and meaningful
absolute tolerances for engineering uncertainty calculations.

Deferred work: optimization, general uncertain arrays, correlated distributions,
and broad NumPy interoperability.

## Development

Use Poetry 2.5.1 and the committed `poetry.lock`:

```sh
poetry sync --with examples
poetry run python -m pytest --cov=worstcase --cov-report=term-missing
poetry run python scripts/check_notebook.py
poetry run python scripts/audit_dependencies.py
poetry build
```

The optional `examples` group provides notebook execution and SciPy without requiring
a notebook web server. Choose your own editor or notebook frontend. CI validates
all five backends, documentation examples, the notebook, and installed wheel/sdist
artifacts containing no unit libraries. Dependabot checks Python dependencies and
GitHub Actions weekly. See [the release checklist](RELEASE.md) for approval gates.
