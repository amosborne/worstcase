# worstcase

For a detailed example of how this software may be leveraged in a true to life example, consider reading this [blog post](https://www.osborneee.com/worstcase/), where the end-to-end measurement uncertainty of a high-side current-sensing circuit is computed.

## What's the worst that could happen?

Professional engineers spend a disproportionate amount of time considering the worst case, especially in fields such as aerospace where the cost of failure can be enormous and therefore the tolerance for technical risk is low.

When delivering hardware to a customer it is typical to also deliver analyses as data products. One such analysis is the worst-case analysis. Hardware performance must be analytically verified to meet requirements for the life of the mission, across all operational environments, with worst-case component variations.

The typical method for performing such an analysis is a spreadsheet like [this one](https://docs.google.com/spreadsheets/d/1OWK2Hds00IrvRUNogDVzHMQhLLowioNIzL4SbS0E3kI/edit#gid=0)... the `worstcase` Python package offers a far more effective solution.

## Usage

At its core, the `worstcase` Python package computes three values: the nominal, the lower bound, and the upper bound. These values may be determind either by Extreme Value, Root-Sum-Square, or Monte Carlo methods.

Input parameters are defined by their range or tolerance, (`param.byrange`, `param.bytol`).

```python
# define the resistor uncertainties
R1 = param.bytol(nom=100 * unit.mohm, tol=0.01, rel=True, tag="R1")
R2 = param.bytol(nom=1.001 * unit.kohm, tol=0.01, rel=True, tag="R2")
R3 = param.bytol(nom=50.5 * unit.kohm, tol=0.01, rel=True, tag="R3")
R4 = param.bytol(nom=1.001 * unit.kohm, tol=0.01, rel=True, tag="R4")
R5 = param.bytol(nom=50.5 * unit.kohm, tol=0.01, rel=True, tag="R5")

# define the amplifier offset voltage
VOS = param.bytol(nom=0 * unit.mV, tol=150 * unit.uV, rel=False, tag="VOS")
```

Derived parameters use a decorator to map worst case input parameters to function arguments (`derive.byev`, `derive.bymc`, or `derive.byrss`).

```python
# define the output voltage
@derive.byev(r1=R1, r2=R2, r3=R3, r4=R4, r5=R5, vos=VOS)
def VO(vbus, iload, r1, r2, r3, r4, r5, vos):
    vp = vbus * r3 / (r2 + r3)
    vn = vp + vos
    vo = vn - (vbus - r1 * iload - vn) * r5 / r4
    return vo

# define the end-to-end uncertainty
@derive.byev(r1=R1, r2=R2, r3=R3, r4=R4, r5=R5, vos=VOS)
def IUNC(r1, r2, r3, r4, r5, vos, vbus, iload):
    vo = VO(vbus, iload, r1, r2, r3, r4, r5, vos)
    return vo / VO(vbus, iload).nom * iload - iload
```

The worst case solution is determined by brute force. If desired, the resulting derived parameter can then be used in the definition of another derived parameter to build up a more complicated analysis.

```python
# calculate at 36V, 1A operating point
VOUT_1A = VO(vbus=36 * unit.V, iload=1 * unit.A, tag="VOUT_1A")
IUNC_1A = IUNC(vbus=36 * unit.V, iload=1 * unit.A, tag="IUNC_1A")

print([VOUT_1A, IUNC_1A])

# [VOUT_1A: 5.045 V (nom), 3.647 V (lb), 6.387 V (ub),
#  IUNC_1A: 0 A (nom), -277 mA (lb), 266 mA (ub)]
```

Parameter units are supported via the default [Pint](https://pypi.org/project/Pint/) `UnitRegistry` object. Results can also be further analyzed for their uncertainty drivers by performing a sensitivity study (`ss()`).

```python
# perform sensitivity study at the 36V, 1A operating point
IUNC_1A_sensitivities = [
    IUNC_1A(tag="IUNC_1A-R1-sens").ss(R1),
    IUNC_1A(tag="IUNC_1A-VOS-sens").ss(VOS),
    IUNC_1A(tag="IUNC_1A-R2-thru-R5-sens").ss([R2, R3, R4, R5]),
]

print(IUNC_1A_sensitivities)

# [IUNC_1A-R1-sens: 0 A (nom), -10 mA (lb), 10 mA (ub),
#  IUNC_1A-VOS-sens: 0 A (nom), -1.53 mA (lb), 1.53 mA (ub),
#  IUNC_1A-R2-thru-R5-sens: 0 A (nom), -265.3 mA (lb), 254.7 mA (ub)]
```

## Installation

`pip install worstcase`
