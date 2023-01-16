import pytest

from worstcase import derive, param, unit


def test_param_byrange():
    p = param.byrange(1.23456, 1, 2, tag="test_byrange", sigfig=3)
    assert str(p) == "test_byrange: 1.23 (nom), 1 (lb), 2 (ub)"
    assert p.nom == 1.23456 and p.lb == 1 and p.ub == 2


def test_param_bytol_absolute():
    p = param.bytol(1, 0.123, False, tag="test_bytol_absolute", sigfig=2)
    assert str(p) == "test_bytol_absolute: 1 (nom), 0.88 (lb), 1.1 (ub)"
    assert p.nom == 1 and p.lb == 0.877 and p.ub == 1.123


def test_param_bytol_relative():
    p = param.bytol(2, 0.16, True, sigfig=2)
    assert str(p) == "2 (nom), 1.7 (lb), 2.3 (ub)"
    assert p.nom == 2 and p.lb == 1.68 and p.ub == 2.32


def test_param_units():
    p = param.bytol(1 * unit.A, 0.1, True)
    assert p.units == unit.A
    assert p.ito(unit("C/s")).nom.m == 1


def test_param_outoforder():
    with pytest.raises(AssertionError):
        param.byrange(0, 1, 1)


def test_derive_byev():
    A = param.byrange(5, 0, 10)
    B = param.bytol(2, 0.1, False)
    C = derive.byev(A, B)(lambda a, b: a + b)
    assert C.nom == 7 and C.lb == 1.9 and C.ub == 12.1
    assert C(a=6).nom == 8 and C(a=6).lb == 7.9 and C(a=6).ub == 8.1
    assert C(a=6, b=2.05) == 8.05


def test_derive_bymc():
    A = param.byrange(5, 0, 10)
    B = param.bytol(2, 0.1, False)
    C = derive.bymc(A, B, n=5000)(lambda a, b: a + b)
    assert (
        C.nom == 7
        and C.lb == pytest.approx(1.9, abs=0.05)
        and C.ub == pytest.approx(12.1, abs=0.05)
    )
    assert (
        C(a=6).nom == 8
        and C(a=6).lb == pytest.approx(7.9, abs=0.05)
        and C(a=6).ub == pytest.approx(8.1, abs=0.05)
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
