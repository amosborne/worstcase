import pytest

from worstcase import derive, param, unit


def test_param_byrange():
    p = param.byrange(1.23456, 1, 2, tag="test_byrange", sigfig=3)
    assert p.nom == 1.23456 and p.lb == 1 and p.ub == 2
    assert str(p) == "test_byrange: 1.23 (nom), 1 (lb), 2 (ub)"
    assert f"{p:L}" == "$1.23 \\,{}^{+0.765}_{-0.235}$"


def test_param_bytol_absolute():
    p = param.bytol(1, 0.123, False, tag="test_bytol_absolute", sigfig=2)
    assert p.nom == 1 and p.lb == 0.877 and p.ub == 1.123
    assert str(p) == "test_bytol_absolute: 1 (nom), 0.88 (lb), 1.1 (ub)"
    assert f"{p:L}" == "$1 \\pm 0.12$"


def test_param_bytol_relative():
    p = param.bytol(2, 0.16, True, sigfig=2)
    assert p.nom == 2 and p.lb == 1.68 and p.ub == 2.32
    assert str(p) == "2 (nom), 1.7 (lb), 2.3 (ub)"
    assert f"{p:L}" == "$2 \\pm 0.32$"


def test_param_units():
    p = param.bytol(1 * unit.A, 0.1, True)
    assert p.units == unit.A
    assert p.ito(unit("C/s")).nom.m == 1


def test_param_pint():
    p = param.bytol(1 * unit.feet / unit.hour, 0.1, True)
    assert p.nom.u == unit.feet / unit.hour
    assert p.nom.u == p.lb.u == p.ub.u
    p.ito(unit.miles / unit.second)
    assert p.nom.u == unit.mile / unit.second
    assert p.nom.u == p.lb.u == p.ub.u
    p.ito_base_units()
    assert p.nom.u == unit.meter / unit.second
    assert p.nom.u == p.lb.u == p.ub.u


def test_param_outoforder():
    with pytest.raises(ValueError):
        param.byrange(0, 1, 1)


def test_param_different_units():
    with pytest.raises(ValueError):
        param.byrange(1 * unit.second, 0.8 * unit.ms, 2 * unit.hour)

<<<<<<< HEAD
=======
    with pytest.raises(ValueError):
        param.byrange(1 * unit.second, 2, 0.5)

>>>>>>> 9463217 (ref: change user-reachable assertions to exceptions)

def test_derive_byev():
    A = param.byrange(5, 0, 10, tag="A")
    B = param.bytol(2, 0.1, False, tag="B")
    C = derive.byev(A, B, tag="C")(lambda a, b: a + b)
    assert C.nom == 7 and C.lb == 1.9 and C.ub == 12.1
    assert C(a=6).nom == 8 and C(a=6).lb == 7.9 and C(a=6).ub == 8.1
    assert C(a=6, b=2.05) == 8.05
    assert C.derivation.nom == {A: A.nom, B: B.nom}
    assert C.derivation.lb == {A: A.lb, B: B.lb}
    assert C.derivation.ub == {A: A.ub, B: B.ub}
    assert str(C) == "C: 7 (nom), 1.9 (lb), 12.1 (ub)"
    assert f"{C:L}" == "$7 \\pm 5.1$"


def test_derive_bymc():
    A = param.byrange(5, 0, 10)
    B = param.bytol(2, 0.1, False)
    C = derive.bymc(A, B, n=5000)(lambda a, b: a + b)
    assert (
        C.nom == 7
        and C.lb == pytest.approx(1.9, abs=0.06)
        and C.ub == pytest.approx(12.1, abs=0.06)
    )
    assert (
        C(a=6).nom == 8
        and C(a=6).lb == pytest.approx(7.9, abs=0.06)
        and C(a=6).ub == pytest.approx(8.1, abs=0.06)
    )
    assert C(a=6, b=2.05) == 8.05


def test_derive_byrss():
    A = param.bytol(1, 2, False)
    B = param.bytol(2, 5, False)
    C = derive.byrss(A, B)(lambda a, b: a * b)
    assert C.nom == 2 and C.ub == pytest.approx(8.40312, abs=1e-5)


def test_derive_byrss_warnasymmetric():
    A = param.byrange(0, 0, 5)
    B = param.bytol(0, 1, False)
    C = derive.byrss(A, B)(lambda a, b: a + b)
    with pytest.warns(UserWarning):
        C()


def test_import_units():
    with pytest.raises(ImportError, match=r"^worstcase no longer*"):
        from worstcase import unit  # noqa F401
