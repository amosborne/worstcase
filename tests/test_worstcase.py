import pytest

import worstcase as wca

wca.config.n = 1000


def test_param():

    # create param by range
    a = wca.param.byrange(1, 0, 3)
    assert a.nom == 1 and a.lb == 0 and a.ub == 3

    # create param by absolute tol
    a = wca.param.bytol(2, 0.1, False)
    assert a.nom == 2 and a.lb == 1.9 and a.ub == 2.1

    # create param by relative tol
    a = wca.param.bytol(2, 0.1, True)
    assert a.nom == 2 and a.lb == 1.8 and a.ub == 2.2

    # create param with out of order bounds
    with pytest.raises(AssertionError):
        wca.param.byrange(0, 1, -1)


def test_derive_onelayer():

    A = wca.param.byrange(5, 0, 10)

    # create a partial binding to extreme value
    @wca.derive.byev(A)
    def Cev(a, b):
        return a * b

    # no argument = build the param (error)
    with pytest.raises(TypeError):
        Cev()

    # rebind, some params = return new param
    Dev = Cev(b=12)
    assert Dev.nom == 60 and Dev.lb == 0 and Dev.ub == 120


def test_evmc_onelayer():

    A = wca.param.byrange(5, 0, 10)
    B = wca.param.bytol(2, 1, False)

    # create a full binding to extreme value
    @wca.derive.byev(A, B)
    def Cev(a, b):
        return a * b

    # create a full binding to monte carlo
    @wca.derive.bymc(A, B)
    def Cmc(a, b):
        return a * b

    assert Cev.nom == 10 and Cmc.nom == 10
    assert Cev.lb == 0 and Cmc.lb == pytest.approx(0, abs=1)
    assert Cev.ub == 30 and Cmc.ub == pytest.approx(30, abs=1)
