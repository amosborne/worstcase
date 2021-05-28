# worstcase

## Overview

`pip install worstcase`

Worst case analysis and sensitivity studies using Extreme Value and/or Monte Carlo methods.

This package coexists alongside far more capable uncertainty analysis and error propagation packages such as [uncertainties](https://pypi.org/project/uncertainties/) (first-order propagation), [soerp](https://pypi.org/project/soerp/) (second-order propagation), and [mcerp](https://pypi.org/project/mcerp/) (Monte Carlo propagation).

This package is designed for engineering applications where worst case analysis computations are often done using the Extreme Value method over single-valued functions while falling back to the Monte Carlo method when the worst case is known to not exist at the extremes. The Extreme Value method is implemented as a brute-force search; the Monte Carlo method is implemented with Latin Hypercube Sampling over a uniform distribution.

This package does not transparently handle parameter covariance across functions. Instead, the graph of parameter dependence must be a tree (and therefore acyclic). Constructing a worst case analysis requires explicit definition of what the underlying parameter dependencies are.

## Usage

```python
import worstcase as wca
import numpy as np

wca.config.n = 2000  # Number of Monte Carlo runs. (default: 5000)
wca.config.sigfig = 4  # Number of significant firues to print. (default: 3)
```

The primitive varying parameter is a *Param*. *Params* are constructed either by tolerance (absolute or relative) or by range. *Params* may be given units from the default unit registry of [Pint](https://pint.readthedocs.io/en/stable/) (default is unitless). A tag for printing is also optional (default is an empty string). A *Param* has the fields *nom* (nominal value), *lb* (lower bound), and *ub* (upper bound).

```python
spd_initial = wca.param.bytol(nom=2, tol=0.1, rel=True, unit=wca.unit("m/s"), tag="v0")
accel = wca.param.byrange(nom=0.2, lb=0.1, ub=0.5, unit=wca.unit("m/s**2"), tag="a")
distance = wca.param.byrange(nom=1, lb=0.8, ub=1.1, unit=wca.unit.km, tag="x")

print([spd_initial, accel, distance])
```
```
[v0: 2 m/s (nom), 1.8 m/s (lb), 2.2 m/s (ub),
 a: 200 mm/s² (nom), 100 mm/s² (lb), 500 mm/s² (ub),
 x: 1 km (nom), 800 m (lb), 1.1 km (ub)]
```

A more complex parameter is built up as a *ParamBuilder*. *ParamBuilders* are constructed by decorating single-valued functions by *ev* (Extreme Value) or *mc* (Monte Carlo). The arguments passed to the *ParamBuilder* are partially bound to the underlying function. A parameter dependency tree can be drawn using the assigned tags; *ParamBuilder* will assume a tag corresponding to the function name as a default.

```python
@wca.param.ev(spd_initial, accel, distance)
def spd_final(v, a, x):
    return np.sqrt(v ** 2 + 2 * a * x)


print(spd_final)
```
```
spd_final (ev)
├── a
├── v0
└── x
```

*ParamBuilder* is a callable and the returned value depends on the arguments supplied. If no arguments are supplied, the parameter is built and a *Param* is returned.

```python
print(spd_final())
```
```
spd_final: 20.1 m/s (nom), 12.78 m/s (lb), 33.24 m/s (ub)
```

Alternatively, the *ParamBuilder* binding to the underlying function can be updated and a new *ParamBuilder* is returned.

```python
spd_final_noaccel = spd_final(a=0 * wca.unit("m/s**2"), tag="spd_noaccel")
print(spd_final_noaccel)
```
```
spd_noaccel (ev)
├── v0
└── x
```

Finally, if the *ParamBuilder* binding is updated such that no arguments are varying parameters then the underlying function will be called to return a single value.

```python
result = spd_final_noaccel(3 * wca.unit("m/s"), x=10 * wca.unit.m)
print(result)
```
```
3.0 meter / second
```

*ParamBuilders* can be used to construct other *ParamBuilders*.

```python
spd_rel = wca.param.bytol(nom=20, tol=1, rel=False, unit=wca.unit("mi/hr"), tag="vrel")


@wca.param.mc(spd_final, spd_rel)
def spd_total(vf, vr):
    return vf + vr


print(spd_total)
print(spd_total())
```
```
spd_total (mc)
├── spd_final (ev)
│   ├── a
│   ├── v0
│   └── x
└── vrel

spd_total: 29.04 m/s (nom), 21.36 m/s (lb), 42.52 m/s (ub)
```

*ParamBuilders* can be modified with the *ss* method to perform a sensitivity study. By supplying a *Param* or *ParamBuilder* (or a list of them), a new *ParamBuilder* is returned where all other varying parameters are set to their nominal value. A few examples below.

```python
accel_sens = spd_total.ss(accel, tag="accel-sens")
print(accel_sens)
print(accel_sens())
```
```
accel-sens (mc)
└── spd_final (ev)
    └── a

accel-sens: 29.04 m/s (nom), 23.23 m/s (lb), 40.62 m/s (ub)
```

```python
accel_distance_sens = spd_total.ss([accel, distance], tag="accel/distance-sens")
print(accel_distance_sens)
print(accel_distance_sens())
```
```
accel/distance-sens (mc)
└── spd_final (ev)
    ├── a
    └── x

accel/distance-sens: 29.04 m/s (nom), 21.75 m/s (lb), 42.16 m/s (ub)
```

```python
finalspd_sens = spd_total.ss(spd_final, tag="finalspd-sens")
print(finalspd_sens)
print(finalspd_sens())
```
```
finalspd-sens (mc)
└── spd_final (ev)
    ├── a
    ├── v0
    └── x

finalspd-sens: 29.04 m/s (nom), 21.73 m/s (lb), 42.18 m/s (ub)
```

```python
relspd_sens = spd_total.ss(spd_rel, tag="relspd-sens")
print(relspd_sens)
print(relspd_sens())
```
```

relspd-sens (mc)
└── vrel

relspd-sens: 29.04 m/s (nom), 28.59 m/s (lb), 29.49 m/s (ub)
```

*Params* implement a few of the Pint *Quantity* methods for usability.

```python
assert spd_total.check("[length]/[time]")
# also try: spd_total.dimensionality

print(spd_total.units)  # meter / second
# also try: spd_total.u

print(spd_total().ito(wca.unit("km/hr")))  # spd_total: 104.5 km/hr (nom)...
# also try: spd_total().ito_base_units()
#           spd_total().ito_reduced_units()
#           spd_total().ito_root_units()
```