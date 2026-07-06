import forallpeople as si
import pytest
import quantities as pq
import unyt
from astropy import units as u
from pint import UnitRegistry

from worstcase import derive, param

unit = UnitRegistry()
si.environment("default", top_level=False)


@pytest.mark.parametrize(
    "amp_unit, volt_unit, ohm_unit",
    [
        pytest.param(1, 1, 1, id="unitless"),
        pytest.param(unit.A, unit.V, unit.ohm, id="pint"),
        pytest.param(si.A, si.V, si.Ohm, id="forallpeople"),
        pytest.param(u.A, u.V, u.Ohm, id="astropy"),
        pytest.param(unyt.A, unyt.V, unyt.Ohm, id="unyt"),
        pytest.param(pq.A, pq.V, pq.Ohm, id="python-quantities"),
    ],
)
class TestUnitizedFunctions:
    """Tests that are run across all supported unit libraries.

    Note that all functions need to have all thiee unites as arguments or
        else pytest will complain. Not all tests use all units.
    """

    def test_param_byrange(self, amp_unit, volt_unit, ohm_unit):
        p = param.byrange(1.23456 * amp_unit, 1 * amp_unit, 2 * amp_unit, sigfig=3)
        assert p.nom == 1.23456 * amp_unit
        assert p.lb == 1 * amp_unit
        assert p.ub == 2 * amp_unit

    def test_param_bytol_absolute(self, amp_unit, volt_unit, ohm_unit):
        p = param.bytol_abs(1 * volt_unit, 0.123 * volt_unit, sigfig=2)
        assert p.nom == 1 * volt_unit
        assert p.lb == 0.877 * volt_unit
        assert p.ub == 1.123 * volt_unit

        p1 = param.bytol(1 * volt_unit, 0.123 * volt_unit, rel=False, sigfig=2)
        assert str(p) == str(p1)

    def test_param_bytol_relative(self, amp_unit, volt_unit, ohm_unit):
        p = param.bytol_rel(2 * ohm_unit, 0.16, sigfig=2)
        assert p.nom == 2 * ohm_unit
        assert p.lb == 1.68 * ohm_unit
        assert p.ub == 2.32 * ohm_unit

        p1 = param.bytol(2 * ohm_unit, 0.16, True, sigfig=2)
        assert str(p) == str(p1)

    def test_param_bytol_relative_invalid_input(
        self, amp_unit, volt_unit, ohm_unit, request
    ):
        if "unitless" in request.node.callspec.id:
            return

        with pytest.raises(ValueError):
            param.bytol_rel(2 * ohm_unit, 0.1 * ohm_unit)

        with pytest.raises(ValueError):
            param.bytol_rel(2 * ohm_unit, 0.1 * amp_unit)

        with pytest.raises(ValueError):
            param.bytol_rel(2, 0.1 * amp_unit)

    def test_param_bytol_absolute_invalid_input(
        self, amp_unit, volt_unit, ohm_unit, request
    ):
        if "unitless" in request.node.callspec.id:
            return

        with pytest.raises(ValueError):
            param.bytol_abs(1 * volt_unit, 0.5)

        with pytest.raises(ValueError):
            param.bytol_abs(1 * volt_unit, 0.5 * amp_unit)

    def test_param_byrange_invalid_input(self, amp_unit, volt_unit, ohm_unit, request):
        if "unitless" in request.node.callspec.id:
            return

        with pytest.raises(ValueError):
            param.byrange(1 * volt_unit, 0.5, 1.5)

        with pytest.raises(ValueError):
            param.byrange(1 * volt_unit, 0.5, 1.5 * volt_unit)

        with pytest.raises(ValueError):
            param.byrange(1 * volt_unit, 0.5 * amp_unit, 1.5 * amp_unit)

    def test_derive_byev(self, amp_unit, volt_unit, ohm_unit):
        A = param.byrange(5 * amp_unit, 0 * amp_unit, 10 * amp_unit, tag="A")
        B = param.bytol(2 * amp_unit, 0.1 * amp_unit, False, tag="B")
        C = derive.byev(A, B, tag="C")(lambda a, b: a + b)
        assert C.nom == 7 * amp_unit
        assert C.lb == 1.9 * amp_unit
        assert C.ub == 12.1 * amp_unit

        assert C(a=6 * amp_unit).nom == 8 * amp_unit
        assert C(a=6 * amp_unit).lb == 7.9 * amp_unit
        assert C(a=6 * amp_unit).ub == 8.1 * amp_unit

        assert C(a=6 * amp_unit, b=2.05 * amp_unit) == 8.05 * amp_unit
        assert C.derivation.nom == {A: A.nom, B: B.nom}
        assert C.derivation.lb == {A: A.lb, B: B.lb}
        assert C.derivation.ub == {A: A.ub, B: B.ub}

    # in case the Monte Carlo sim is "unlucky" w/o having to crank up iterations
    @pytest.mark.flaky(reruns=3, only_rerun=["AssertionError"])
    def test_derive_bymc(self, amp_unit, volt_unit, ohm_unit):
        A = param.byrange(5 * ohm_unit, 0 * ohm_unit, 10 * ohm_unit)
        B = param.bytol(2 * ohm_unit, 0.1 * ohm_unit, False)
        C = derive.bymc(A, B, n=5000)(lambda a, b: a + b)
        assert (
            C.nom == 7 * ohm_unit
            and C.get_val(C.lb) == pytest.approx(1.9, abs=0.06)
            and C.get_val(C.ub) == pytest.approx(12.1, abs=0.06)
        )
        assert (
            C(a=6 * ohm_unit).nom == 8 * ohm_unit
            and C.get_val(C(a=6 * ohm_unit).lb) == pytest.approx(7.9, abs=0.06)
            and C.get_val(C(a=6 * ohm_unit).ub) == pytest.approx(8.1, abs=0.06)
        )
        assert C(a=6 * ohm_unit, b=2.05 * ohm_unit) == 8.05 * ohm_unit

    def test_derive_byrss(self, amp_unit, volt_unit, ohm_unit):
        A = param.bytol(1 * amp_unit, 2 * amp_unit, False)
        B = param.bytol(2 * volt_unit, 5 * volt_unit, False)
        C = derive.byrss(A, B)(lambda a, b: a * b)
        assert C.nom == 2 * amp_unit * volt_unit
        assert C.get_val(C.ub) == pytest.approx(8.40312, abs=1e-5)

    def test_complex_byev(self, amp_unit, volt_unit, ohm_unit):
        # define the resistor uncertainties
        R1 = param.bytol(nom=100e-3 * ohm_unit, tol=0.01, rel=True, tag="R1")
        R2 = param.bytol(nom=1.001e3 * ohm_unit, tol=0.01, rel=True, tag="R2")
        R3 = param.bytol(nom=50.5e3 * ohm_unit, tol=0.01, rel=True, tag="R3")
        R4 = param.bytol(nom=1.001e3 * ohm_unit, tol=0.01, rel=True, tag="R4")
        R5 = param.bytol(nom=50.5e3 * ohm_unit, tol=0.01, rel=True, tag="R5")

        # define the amplifier offset voltage
        VOS = param.bytol(
            nom=0 * volt_unit, tol=150e-6 * volt_unit, rel=False, tag="VOS"
        )

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

        # calculate at 36V, 1A operating point
        VOUT_1A = VO(vbus=36 * volt_unit, iload=1 * amp_unit, tag="VOUT_1A")
        IUNC_1A = IUNC(vbus=36 * volt_unit, iload=1 * amp_unit, tag="IUNC_1A")

        get_val = VOUT_1A.get_val
        get_units = VOUT_1A.get_units

        assert get_val(VOUT_1A.nom) == pytest.approx(5.045, abs=1e-3)
        assert get_val(VOUT_1A.lb) == pytest.approx(3.647, abs=1e-3)
        assert get_val(VOUT_1A.ub) == pytest.approx(6.387, abs=1e-3)
        assert get_units(VOUT_1A.ub) == volt_unit
        assert IUNC_1A.nom == 0 * amp_unit
        assert get_val(IUNC_1A.lb) == pytest.approx(-0.277, abs=1e-3)
        assert get_val(IUNC_1A.ub) == pytest.approx(0.266, abs=1e-3)
        assert get_units(IUNC_1A.ub) == amp_unit

    def test_param_equivalency(self, amp_unit, volt_unit, ohm_unit, request):
        A = param.byrange(5 * amp_unit, 0 * amp_unit, 10 * amp_unit, tag="A")
        B = param.bytol(2 * amp_unit, 0.1 * amp_unit, False, tag="B")
        C = derive.byev(A, B, tag="C")(lambda a, b: a + b)
        assert A.equivalent(
            param.byrange(5 * amp_unit, 0 * amp_unit, 10 * amp_unit, tag="A")
        )
        assert A.equivalent(param.byrange(5 * amp_unit, 0 * amp_unit, 10 * amp_unit))
        assert not A.equivalent(
            param.byrange(6 * amp_unit, 0 * amp_unit, 10 * amp_unit)
        )
        assert not A.equivalent(
            param.byrange(5 * amp_unit, 1 * amp_unit, 10 * amp_unit)
        )
        assert not A.equivalent(param.byrange(5 * amp_unit, 0 * amp_unit, 9 * amp_unit))
        assert B.equivalent(
            param.byrange(2 * amp_unit, 1.9 * amp_unit, 2.1 * amp_unit, tag="B")
        )
        assert B.equivalent(param.byrange(2 * amp_unit, 1.9 * amp_unit, 2.1 * amp_unit))
        assert C.equivalent(
            param.byrange(7 * amp_unit, 1.9 * amp_unit, 12.1 * amp_unit, tag="C")
        )
        assert C.equivalent(
            param.byrange(7 * amp_unit, 1.9 * amp_unit, 12.1 * amp_unit)
        )

        if "unitless" in request.node.callspec.id:
            return

        assert not B.equivalent(
            param.byrange(2 * volt_unit, 1.9 * volt_unit, 2.1 * volt_unit)
        )
        assert not A.equivalent(
            param.byrange(5 * volt_unit, 0 * volt_unit, 10 * volt_unit)
        )


def test_param_str_repr():
    p1 = param.bytol(2, 0.16, True, sigfig=2)
    p2 = param.bytol(1, 0.123, False, tag="absolute", sigfig=2)
    p3 = param.byrange(1.23456, 1, 2, tag="range", sigfig=3)

    assert str(p1) == "2 (nom), 1.7 (lb), 2.3 (ub)"
    assert str(p2) == "absolute: 1 (nom), 0.88 (lb), 1.1 (ub)"
    assert str(p3) == "range: 1.23 (nom), 1 (lb), 2 (ub)"


def test_param_latex_repr():
    p1 = param.bytol(2, 0.16, True, sigfig=2)
    p2 = param.bytol(1, 0.123, False, tag="absolute", sigfig=2)
    p3 = param.byrange(1.23456, 1, 2, tag="range", sigfig=3)

    assert f"{p1:L}" == "$2 \\pm 0.32$"
    assert f"{p2:L}" == "$1 \\pm 0.12$"
    assert f"{p3:L}" == "$1.23 \\,{}^{+0.765}_{-0.235}$"


@pytest.mark.parametrize(
    "base_unit, scaled_unit",
    [
        pytest.param(1, 1e-3, id="unitless"),
        pytest.param(unit.V, unit.mV, id="pint"),
        pytest.param(u.V, u.mV, id="astropy"),
        pytest.param(unyt.V, unyt.mV, id="unyt"),
        pytest.param(pq.V, pq.mV, id="python-quantities"),
    ],
)
def test_param_similar_units(base_unit, scaled_unit):
    param.byrange(4 * base_unit, 1 * scaled_unit, 5 * base_unit)
    param.bytol_abs(4 * base_unit, 1 * scaled_unit)

    with pytest.raises(ValueError):
        param.byrange(4000 * scaled_unit, 5 * base_unit, 6000 * scaled_unit)


def test_param_outoforder():
    with pytest.raises(ValueError):
        param.byrange(0, 1, 1)


def test_equality_different_units():
    A = param.bytol(2 * u.A, 0.1, rel=True)
    B = param.bytol(2 * si.A, 0.1, rel=True)
    C = param.bytol(2, 0.1, rel=True)
    assert not A.equivalent(B)
    assert not B.equivalent(A)
    assert not C.equivalent(A)
    assert not C.equivalent(B)


def test_derive_byrss_warnasymmetric():
    A = param.byrange(0, 0, 5)
    B = param.bytol(0, 1, False)
    C = derive.byrss(A, B)(lambda a, b: a + b)
    with pytest.warns(UserWarning):
        C()


def test_import_units():
    with pytest.raises(ImportError, match=r"^worstcase no longer"):
        from worstcase import unit  # noqa F401
