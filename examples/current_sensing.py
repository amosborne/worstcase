"""Current-sensing example shared by the notebook and regression tests."""

import numpy as np
from pint import UnitRegistry

from worstcase import derive, param

unit = UnitRegistry()
R1 = param.bytol_rel(100 * unit.mohm, 0.01)
R2 = param.bytol_rel(1.001 * unit.kohm, 0.01)
R3 = param.bytol_rel(50.5 * unit.kohm, 0.01)
R4 = param.bytol_rel(1.001 * unit.kohm, 0.01)
R5 = param.bytol_rel(50.5 * unit.kohm, 0.01)
VOS = param.bytol_abs(0 * unit.mV, 150 * unit.uV)
inputs = (R1, R2, R3, R4, R5, VOS)


def output_voltage(vbus, iload, r1, r2, r3, r4, r5, vos):
    """Works with expressions or concrete quantities."""
    vp = vbus * r3 / (r2 + r3)
    vn = vp + vos
    return vn - (vbus - r1 * iload - vn) * r5 / r4


def operating_point(vbus=36 * unit.V, iload=1 * unit.A):
    """Return output-voltage and current-error expressions."""
    bus, load = param.constant(vbus), param.constant(iload)
    voltage = output_voltage(bus, load, *inputs)
    current_error = (voltage / voltage.nom - 1) * load
    return voltage, current_error


@derive
def common_mode_gain(r1, r2, r3, r4, r5, vos):
    """A scalar result from an internal array-based linear solve."""
    from scipy.linalg import solve

    rows, outputs = [], []
    for bus, load in [(0, 0), (1, 10), (3, 20)]:
        bus, load = bus * unit.V, load * unit.A
        vo = output_voltage(bus, load, r1, r2, r3, r4, r5, vos)
        vdm, vcm = r1 * load, bus - r1 * load / 2
        outputs.append(vo.to("V").magnitude)
        rows.append([vcm.to("V").magnitude, vdm.to("V").magnitude, 1])
    return float(solve(np.array(rows), np.array(outputs))[0])


if __name__ == "__main__":
    voltage, current_error = operating_point()
    print(format(voltage.derive(), ".4g"))
    print(format(current_error.derive(), ".4g"))
    print(format(common_mode_gain(*inputs).derive(), ".4g"))
