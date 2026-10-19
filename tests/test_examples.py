import re
from pathlib import Path

import pytest

from examples.current_sensing import (
    R1,
    R2,
    R3,
    R4,
    R5,
    VOS,
    common_mode_gain,
    inputs,
    operating_point,
    unit,
)


def test_current_sensing_results_and_sensitivities():
    voltage, current_error = operating_point()
    vo, error = voltage.derive(), current_error.derive()
    assert vo.nom.to("V").magnitude == pytest.approx(5.045, abs=1e-3)
    assert vo.lb.to("V").magnitude == pytest.approx(3.647, abs=1e-3)
    assert vo.ub.to("V").magnitude == pytest.approx(6.387, abs=1e-3)
    assert error.lb.to("A").magnitude == pytest.approx(-0.277, abs=1e-3)
    assert error.ub.to("A").magnitude == pytest.approx(0.266, abs=1e-3)
    resistor = current_error.ss(R1)
    offset = current_error.ss(VOS)
    network = current_error.ss([R2, R3, R4, R5])
    assert resistor.lb.to("mA").magnitude == pytest.approx(-10, abs=0.01)
    assert offset.ub.to("mA").magnitude == pytest.approx(1.53, abs=0.01)
    assert network.lb.to("mA").magnitude == pytest.approx(-265.3, abs=0.1)
    assert network.ub.to("mA").magnitude == pytest.approx(254.7, abs=0.1)
    low_current = operating_point(iload=0.05 * unit.A)[1].derive()
    assert low_current.nom.to("A").magnitude == pytest.approx(0, abs=1e-12)


def test_internal_linear_solve_matches_resistor_network():
    gain = common_mode_gain(*inputs)
    result = gain.derive()
    # Analytical common-mode coefficient from the voltage equation.
    ratio = R3.nom / (R2.nom + R3.nom)
    expected = (ratio * (1 + R5.nom / R4.nom) - R5.nom / R4.nom).magnitude
    assert result.nom == pytest.approx(expected, abs=1e-12)
    assert result.lb < result.nom < result.ub
    assert gain.ss([R2, R3, R4, R5]).lb == pytest.approx(result.lb)


def test_readme_python_examples_execute():
    namespace = {}
    readme = Path(__file__).resolve().parents[1] / "README.md"
    for source in re.findall(r"```python\n(.*?)```", readme.read_text(), re.DOTALL):
        exec(compile(source, str(readme), "exec"), namespace)
