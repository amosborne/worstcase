# worstcase

For a worked engineering example, see my
[blog post](https://www.osborneee.com/worstcase/) on the end-to-end measurement
uncertainty of a high-side current-sensing circuit, and the
[updated notebook](https://github.com/amosborne/worstcase/blob/release/1.0/example.ipynb).

## What's the worst that could happen?

Professional engineers spend a disproportionate amount of time considering the worst
case, especially in fields such as aerospace where the cost of failure can be
enormous and therefore the tolerance for technical risk is low.

When delivering hardware to a customer it is typical to also deliver analyses as
data products. One such analysis is the worst-case analysis. Hardware performance
must be analytically verified to meet requirements for the life of the mission,
across all operational environments, with worst-case component variations.

The typical method for performing such an analysis is a spreadsheet like
[this one](https://docs.google.com/spreadsheets/d/1OWK2Hds00IrvRUNogDVzHMQhLLowioNIzL4SbS0E3kI/edit#gid=0).
The `worstcase` Python package lets you express these calculations as reusable
formulas, with or without units. Write a formula once, then evaluate it using
corner analysis, Monte Carlo sampling, or root-sum-square analysis. Repeated uses
of the same input stay linked throughout the calculation.

The reported range depends on the analysis method; it is **not a general guarantee
of worst-case bounds**. The examples below show both the capabilities and the limits.

## Install

Requires Python 3.11 or later; tested on Python 3.11–3.14, Windows and Linux.

```sh
pip install worstcase
```

Upgrading from 0.x? See the
[migration notes below](#upgrading-from-0x).

Plain numbers need no unit package. For the Pint examples below, also install:

```sh
pip install pint
```

Pint, forallpeople, Astropy, unyt, and python-quantities are supported. Install the
library you want separately; there are no unit-library extras to select.

## A first calculation

A 5 V signal with 2% tolerance feeds an amplifier with gain 2 ± 0.1. Using plain
numbers, with volts understood:

```python
from worstcase import param

signal = param.bytol_rel(5, 0.02)  # 0.02 means 2%, not 0.02%.
gain = param.bytol_abs(2, 0.1)
output = signal * gain

result = output.derive(tag="output")  # Default: corner analysis ("ev").
print(f"{result:.4g}")
print(f"nominal={result.nom:.2f}, lower={result.lb:.2f}, upper={result.ub:.2f}")
```

```text
output: 10 (nom), 9.31 (lb), 10.71 (ub)
nominal=10.00, lower=9.31, upper=10.71
```

Build a formula with ordinary arithmetic, then call `.derive()` to calculate its
bounds. `.nom` on a formula calculates only its nominal value. An evaluated result
stores `.nom`, `.lb`, and `.ub`; accessing or displaying them does not rerun the model.

The Python examples below build on earlier imports and definitions. Each `text`
block shows the preceding example's printed output.

## Describe your inputs

Choose a relative tolerance, an absolute tolerance, explicit asymmetric bounds,
or a fixed value:

```python
relative = param.bytol_rel(100, 0.05, tag="relative")
absolute = param.bytol_abs(100, 2, tag="absolute")
asymmetric = param.byrange(100, 97, 106, tag="asymmetric")
fixed = param.constant(100, tag="fixed")

for value in (relative, absolute, asymmetric, fixed):
    print(f"{value:.0f}")
```

```text
relative: 100 (nom), 95 (lb), 105 (ub)
absolute: 100 (nom), 98 (lb), 102 (ub)
asymmetric: 100 (nom), 97 (lb), 106 (ub)
fixed: 100 (nom), 100 (lb), 100 (ub)
```

Bounds must satisfy `lower <= nominal <= upper`. Tolerances must be nonnegative;
relative tolerances use the magnitude of the nominal value, so a negative nominal
is supported. Use an absolute tolerance or explicit bounds for a nominal of zero.
With units, absolute tolerances and bounds must have compatible dimensions;
relative tolerances must be dimensionless, such as `0.01` or Pint's `1 * u.percent`.

### Reuse an input when it is the same physical quantity

```python
x = param.bytol_abs(10, 1)
y = param.bytol_abs(10, 1)  # A separate input with the same specification.

print(f"{(x - x).derive(tag='same input'):.1f}")
print(f"{(x - y).derive(tag='separate inputs'):.1f}")
print(f"{(x / x).derive(tag='same ratio'):.1f}")
```

```text
same input: 0.0 (nom), 0.0 (lb), 0.0 (ub)
separate inputs: 0.0 (nom), -2.0 (lb), 2.0 (ub)
same ratio: 1.0 (nom), 1.0 (lb), 1.0 (ub)
```

Separate parameters vary independently. Reuse the original parameter or formula
wherever that same quantity appears. Arbitrary correlations between separate
inputs are not supported. Keep intermediate **formulas** in your calculation;
evaluated results cannot participate in formula arithmetic. Creating a new
parameter from a result's bounds discards its original input dependencies.

## Work with units: a voltage divider

Two 10 kΩ resistors, each with 1% tolerance, divide a fixed 5 V supply:

```python
from pint import UnitRegistry

u = UnitRegistry()
r1 = param.bytol_rel(10 * u.kohm, 0.01)
r2 = param.bytol_rel(10 * u.kohm, 0.01)
vin = param.constant(5 * u.V)
vout = vin * r2 / (r1 + r2)
divider_result = vout.derive(tag="vout")

print(f"{divider_result:.4g}")
print(f"lower in mV: {divider_result.lb.to('mV').magnitude:.1f}")
```

```text
vout: 2.5 V (nom), 2.475 V (lb), 2.525 V (ub)
lower in mV: 2475.0
```

Result values are native quantities: use your unit library's conversion methods,
such as Pint's `.to(...)`. Bounds are converted to the nominal result's units.
Use one unit library and one registry per calculation, and configure it before
creating inputs. Unit labels and formatting follow the chosen library.

Formulas accept scalar constants automatically: `r1 * (5 * u.V)` works. When a
quantity comes first, use `param.constant(5 * u.V) * r1` for portable behavior;
some libraries reject `quantity * uncertain_input` before worstcase can handle it.

Other supported libraries work with the same `param` API. For example, with
Astropy installed (`pip install astropy`):

```python
from astropy import units as au

distance = param.byrange(10 * au.m, 950 * au.cm, 10.5 * au.m)
duration = param.constant(2 * au.s)
speed = (distance / duration).derive(tag="speed")
print(f"{speed:.2f}")
```

```text
speed: 5.00 m / s (nom), 4.75 m / s (lb), 5.25 m / s (ub)
```

For the remaining libraries, construct inputs with their native quantities, such
as `10 * unyt.V`, `10 * quantities.V`, or `10 * forallpeople.V` after configuring
forallpeople's environment. Plain numbers can scale quantities, but adding a
bare number to a dimensional quantity does not supply the missing units.

## Choose an analysis method

Here are three analyses of the same sum of two independent inputs:

```python
a = param.bytol_abs(10, 1)
b = param.bytol_abs(20, 2)
total = a + b

ev = total.derive(method="ev", tag="EV")
mc = total.derive(method="mc", n=1000, seed=42, tag="MC")
rss = total.derive(method="rss", tag="RSS")
for analysis in (ev, mc, rss):
    print(f"{analysis:.2f}; evaluations={analysis.n_evaluations}")
```

```text
EV: 30.00 (nom), 27.00 (lb), 33.00 (ub); evaluations=5
MC: 30.00 (nom), 27.04 (lb), 32.95 (ub); evaluations=1001
RSS: 30.00 (nom), 27.76 (lb), 32.24 (ub); evaluations=5
```

| Method | What it evaluates | How to interpret the result |
|---|---|---|
| `ev` (default) | Every combination of input endpoints, plus nominal | Corner envelope. Useful when extrema occur at corners, for example a formula monotonic in each input over the whole range. Interior extrema can be missed. |
| `mc` | `n` Latin hypercube samples across independent uniform input ranges, plus nominal | Observed sample envelope. It is neither a guaranteed bound nor a confidence interval. |
| `rss` | Both endpoints of one input at a time, with other inputs at nominal; combines deviations in quadrature | Endpoint sensitivity approximation for independent contributions. Input tolerances are not automatically standard deviations. |

With `k` varying inputs, EV needs `2**k + 1` model evaluations and RSS needs
`2*k + 1`. MC needs `n + 1`. Fixed inputs do not increase these counts; a fully
fixed calculation evaluates once. EV becomes expensive quickly as inputs grow.
`n` must be positive; `n` and `seed` apply only to MC. Set a seed for repeatability
with the same formula and dependency versions. Exact sampled extrema may change
across dependency versions.

### Corner analysis can miss a true minimum

```python
position = param.byrange(1, 0, 2)
bowl = (position - 0.5) ** 2
print(f"{bowl.derive(tag='corner result'):.2f}")
print("Actual minimum: 0.00 at position=0.5")
```

```text
corner result: 0.25 (nom), 0.25 (lb), 2.25 (ub)
Actual minimum: 0.00 at position=0.5
```

EV checks positions 0, 1 (nominal), and 2. It never checks 0.5. MC may get closer
but still does not prove a global minimum. Neither method establishes that a model
is valid between evaluated points; an unsampled singularity can also go undetected.

### RSS can miss interactions entirely

```python
p = param.byrange(0, -1, 1)
q = param.byrange(0, -1, 1)
interaction = p * q
print(f"{interaction.derive(method='ev', tag='EV'):.1f}")
print(f"{interaction.derive(method='rss', tag='RSS'):.1f}")
```

```text
EV: 0.0 (nom), -1.0 (lb), 1.0 (ub)
RSS: 0.0 (nom), 0.0 (lb), 0.0 (ub)
```

Changing either input alone leaves the product at zero; changing both does not.
RSS therefore reports no variation here. It warns about asymmetric endpoint
responses, but **no warning does not mean a formula is linear or RSS is suitable**.
For a linear sum, the RSS example above combines deviations of 1 and 2 as
`sqrt(1**2 + 2**2)`; no probability or confidence level follows from that operation.

## Isolate a tolerance contribution

Use `.ss(...)` to vary selected inputs and hold everything else at nominal:

```python
print(f"{vout.ss(r1, tag='R1 only'):.4f}")
print(f"{vout.ss([r1, r2], tag='both resistors'):.4f}")
print(f"{vout.ss([], tag='all fixed'):.4f}")
```

```text
R1 only: 2.5000 V (nom), 2.4876 V (lb), 2.5126 V (ub)
both resistors: 2.5000 V (nom), 2.4750 V (lb), 2.5250 V (ub)
all fixed: 2.5000 V (nom), 2.5000 V (lb), 2.5000 V (ub)
```

Select a single parameter, a list, or an intermediate formula to vary all the
inputs used by that formula. Selections must belong to the calculation.
Sensitivity studies also accept `method="mc", n=..., seed=...` or `method="rss"`.
Individual contributions need not add up to the full range of a nonlinear model.

## Inspect the inputs behind a bound

EV and MC record the actual input assignment that produced each reported bound:

```python
upper_inputs = divider_result.derivation.ub
print(f"R1 at upper output: {upper_inputs[r1].to('ohm').magnitude:.0f} ohm")
print(f"R2 at upper output: {upper_inputs[r2].to('ohm').magnitude:.0f} ohm")
print(f"upper output: {divider_result.ub.to('V').magnitude:.3f} V")
print(f"RSS upper assignment: {rss.derivation.ub}")
```

```text
R1 at upper output: 9900 ohm
R2 at upper output: 10100 ohm
upper output: 2.525 V
RSS upper assignment: None
```

Assignments are read-only maps keyed by the **original parameter objects**.
`.derivation.nom` records nominal inputs; `.derivation.lb` and `.derivation.ub`
record bound inputs, including fixed values. A tie keeps the first assignment
encountered, with nominal considered first. RSS combines separate evaluations,
so its bound assignments are `None`.

## Use your own functions and numerical calculations

Decorate a function with `@derive` when ordinary formula arithmetic is not enough.
The function receives concrete values when the calculation runs. It can use
conditionals, arrays, or a numerical solver internally, but must return one finite
real number or supported scalar quantity.

```python
from worstcase import derive

@derive
def clipped_amplifier(value, gain=2, *, ceiling=5):
    return min(value * gain, ceiling)

input_signal = param.byrange(2, 1, 3)
clipped = clipped_amplifier(input_signal, ceiling=5)
print(f"{clipped.derive(tag='clipped'):.1f}")
print(f"ordinary call: {clipped_amplifier(2):.1f}")
```

```text
clipped: 4.0 (nom), 2.0 (lb), 5.0 (ub)
ordinary call: 4.0
```

Calls with ordinary values remain ordinary calls. Defaults, positional arguments,
and keyword arguments work normally. Lists, tuples, and dictionary values can
contain uncertain inputs. Here a function calculates parallel resistance using
an internal NumPy array, with all resistance values expressed in ohms:

```python
import numpy as np

@derive
def parallel(resistances):
    values = np.asarray(resistances)
    return 1 / np.sum(1 / values)

resistors = [param.bytol_rel(100, 0.1), param.bytol_rel(100, 0.1)]
equivalent = parallel(resistors)
print(f"{equivalent.derive(tag='parallel ohms'):.1f}")
```

```text
parallel ohms: 50.0 (nom), 45.0 (lb), 55.0 (ub)
```

Keep functions deterministic and free of external side effects. They may run many
times, once per evaluated input assignment. A failed evaluation raises an
exception; samples are not silently discarded. If using SciPy or another optional
library inside a function, install it separately and check its solver converged.

### Supported arithmetic and NumPy functions

Formulas support `+`, `-`, `*`, `/`, `**`, unary `+`/`-`, and `abs(...)`.
Supported scalar NumPy functions are `sqrt`, `square`, `exp`, `log`, `log10`,
`sin`, `cos`, and `tan`, plus the corresponding arithmetic ufuncs:

```python
decay_time = param.byrange(2, 1, 3)
remaining = np.exp(-decay_time)
print(f"{remaining.derive(tag='remaining fraction'):.4g}")
```

```text
remaining fraction: 0.1353 (nom), 0.04979 (lb), 0.3679 (ub)
```

Arguments must satisfy the mathematical domain and the unit library's rules.
General NumPy array operations, reductions on uncertain arrays, and Python's
`math` functions are not supported directly on uncertain formulas. Wrap a
scalar-output calculation with `@derive` as shown above. Likewise, put a
conditional inside a wrapped function instead of writing `if uncertain_input:`.

## Labels, formatting, and notebooks

Labels are inferred from assignments when source code is available. Use `tag=...`
for a dependable label, or `tag=""` to omit one. A label affects display only;
matching labels never link separate inputs. If inference is unavailable, the
label is empty.

```python
reading = param.byrange(10, 9, 12, tag="sensor")
print(f"{reading:.1f}")
print(f"{reading.derive(tag=''):.2e}")
print(f"{reading:.3gL}")
```

```text
sensor: 10.0 (nom), 9.0 (lb), 12.0 (ub)
1.00e+01 (nom), 9.00e+00 (lb), 1.20e+01 (ub)
$10\,{}^{+2}_{-1}$
```

Standard numeric format specifications such as `.3g`, `.2f`, and `.2e` apply to
each value. The default is `.4g`; append `L` for LaTeX. Formatting rounds the
display only. Notebook displays of parameters and evaluated results use LaTeX
automatically; showing an unevaluated formula does not calculate its bounds.

## Capabilities and limits at a glance

| Supported | Boundary to keep in mind |
|---|---|
| Finite real scalar inputs and results | No uncertain arrays, complex values, booleans, NaN, or infinity |
| Repeated references to the same uncertain input | No arbitrary correlation matrices or correlated distributions |
| Relative, absolute, and asymmetric input ranges | No built-in Gaussian, lognormal, or other configurable probability distributions |
| Five optional unit libraries | One library/registry per calculation; use linear units, rather than affine temperature or logarithmic uncertainty models |
| Arithmetic, selected NumPy functions, and wrapped functions | No automatic global optimization, interval proof, or general NumPy integration |
| Bound input assignments for EV and MC | Recorded extrema reflect only evaluated points; RSS bounds have no joint assignment |
| Stable evaluated results | Retain formulas to preserve dependencies when building larger calculations |

Parameter specifications cannot be edited in place. Create new inputs and rebuild
the formula to study a different specification. Input quantities are copied when
captured; avoid changing shared unit registries or global state during an analysis.
All evaluated outputs must remain finite and dimensionally compatible with the
nominal result. Numerical precision and the validity of your engineering model
still matter, regardless of the method chosen.

For a larger worked circuit, see the
[current-sensing notebook](https://github.com/amosborne/worstcase/blob/release/1.0/example.ipynb).
Report problems or request features through
[GitHub issues](https://github.com/amosborne/worstcase/issues).

## Upgrading from 0.x

1.0 is a breaking release. Choose an analysis method when evaluating the final
formula, and keep intermediate formulas to preserve shared input dependencies.

| Before | Now |
|---|---|
| `from worstcase import unit` | Create units with your preferred library, such as Pint's `UnitRegistry()` |
| `param.bytol(nom, tol, rel=True)` | `param.bytol_rel(nom, tol)` |
| `param.bytol(nom, tol, rel=False)` | `param.bytol_abs(nom, tol)` |
| `@derive.byev(...)` | Ordinary formula arithmetic, or `@derive` with explicit input arguments |
| `@derive.bymc(...)` / `@derive.byrss(...)` | `formula.derive(method="mc", n=1000, seed=42)` / `formula.derive(method="rss")` |
| `C()` / `C.lb` | `result = formula.derive(); result.lb` |
| `sigfig=3` | `format(result, ".3g")` |
| `C(tag="study").ss(A)` | `formula.ss(A, tag="study")` |
| Pint-specific `ito` / `apply` methods | Convert concrete quantities in your code or inside a wrapped function |

Result bound assignments are keyed by the original input objects; RSS has no
joint lower/upper assignment. MC records the sampled extrema, including nominal.

## Contributing

Development uses Poetry 2.5.1 and the committed lockfile. From a checkout:

```sh
poetry sync --with examples
poetry run python -m pytest
poetry run python scripts/check_notebook.py
poetry build
```
