from types import SimpleNamespace

import pytest


@pytest.fixture(
    params=["plain", "pint", "forallpeople", "astropy", "unyt", "quantities"]
)
def backend(request):
    name = request.param
    if name == "plain":
        v = a = 1
        mv, percent = 0.001, 0.01
    elif name == "pint":
        from pint import UnitRegistry

        u = UnitRegistry()
        v, a, mv, percent = u.V, u.A, u.mV, u.percent
    elif name == "forallpeople":
        import forallpeople as u

        u.environment("default", top_level=False)
        v, a, mv, percent = u.V, u.A, 0.001 * u.V, 0.01
    else:
        if name == "astropy":
            from astropy import units as u
        elif name == "unyt":
            import unyt as u
        else:
            import quantities as u
        v, a, mv, percent = u.V, u.A, u.mV, u.percent
    return SimpleNamespace(name=name, V=v, A=a, mV=mv, percent=percent)
