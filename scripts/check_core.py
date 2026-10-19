"""Run in an installation containing only runtime dependencies."""

import sys
from importlib.util import find_spec

from worstcase import derive, param

for library in ("pint", "forallpeople", "astropy", "unyt", "quantities"):
    assert find_spec(library) is None, library
    assert library not in sys.modules

x = param.byrange(2, 1, 3)
assert (x - x).derive().ub == 0
assert (x + 1).derive().lb == 2
assert (x + 1).derive(method="mc", n=10, seed=42).n_evaluations == 11


@derive
def square(value):
    return value**2


assert square(x).derive().ub == 9
print("Core-only installation passed")
