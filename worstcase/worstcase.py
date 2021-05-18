from functools import partial
from itertools import product
from random import sample

import numpy as np
from matplotlib import pyplot as plt  # type: ignore
from scipy.stats import uniform  # type: ignore


class Param:
    def __init__(self, arg, nom, lb, ub, **kwargs):
        """
        A randomly varying parameter with a nominal and extreme values.
            arg:   variable name, passed as **kwargs[arg] to functions by mc/ev
            nom:   nominal value, computed as the mean of an mc distribution
            lb:    lower bound, computed as the min of an mc distribution
            ub:    upper bound, computed as the max of an mc distribution
            dist:  random distribution factory function g(param),
                   returns a random distribution function f(n),
                   returns n random values (optional, default uniform)
            name:  parameter description (optional, default "")
            form:  format string (optional, default "{}")
            scale: format scale factor (optional, default 1)
        """
        assert lb <= nom <= ub

        self.arg = arg
        self.nom = nom
        self.lb = lb
        self.ub = ub

        if "dist" in kwargs:
            self.dist = kwargs["dist"](self)
        else:
            self.dist = lambda n: uniform.rvs(loc=lb, scale=ub - lb, size=n)

        self.name = kwargs.get("name", "")
        self.form = kwargs.get("form", "{}")
        self.scale = kwargs.get("scale", 1)

    @staticmethod
    def byrange(arg, nom, lb, ub, **kwargs):
        """Static constructor by range."""
        return Param(arg, nom, lb, ub, **kwargs)

    @staticmethod
    def bytol(arg, nom, tol, rel, **kwargs):
        """Static constructor by relative or absolute tolerance."""
        tol = nom * tol if rel else tol
        return Param(arg, nom, nom - tol, nom + tol, **kwargs)

    def plot(self, n, bins, ax=None):
        """Plot a normalized histogram of distribution."""
        ax = plt.gca() if ax is None else ax
        ax.hist(
            self.dist(n), bins=bins, density=True, histtype="stepfilled", label=self.arg
        )

        ylim = ax.get_ylim()
        ax.plot([self.lb] * 2, ylim, "r-.", label=None)
        ax.plot([self.nom] * 2, ylim, "r-.", label=None)
        ax.plot([self.ub] * 2, ylim, "r-.", label=None)

    def __repr__(self):
        """Print string representation."""
        name = "{} ({}): ".format(self.name, self.arg)
        lb = self.form.format(self.lb * self.scale)
        nom = self.form.format(self.nom * self.scale)
        ub = self.form.format(self.ub * self.scale)
        return name + lb + " < " + nom + " < " + ub


def mc(params, arg, n, **p_kwargs):
    """Monte Carlo analysis."""

    def decorator(func):
        def wrapper(*args, ss=None, **f_kwargs):
            if ss is not None:
                if not isinstance(ss, list):
                    ss = [ss]
                nom_args = {p.arg: p.nom for p in params if p not in ss}
            else:
                nom_args = {}
                ss = params

            param_args = {p.arg: p.dist(n) for p in ss}
            param_args = [dict(zip(param_args, t)) for t in zip(*param_args.values())]
            if param_args:
                mc_result = [
                    func(*args, **f_kwargs, **p, **nom_args) for p in param_args
                ]
            else:
                mc_result = [func(*args, **f_kwargs, **nom_args)] * n

            def mc_dist(param):
                return partial(sample, mc_result)

            return Param(
                arg=arg,
                nom=np.mean(mc_result),
                lb=min(mc_result),
                ub=max(mc_result),
                dist=mc_dist,
                **p_kwargs
            )

        return wrapper

    return decorator


def ev(params, arg, **p_kwargs):
    """Extreme Value analysis."""

    def decorator(func):
        def wrapper(*args, ss=None, **f_kwargs):
            if ss is not None:
                if not isinstance(ss, list):
                    ss = [ss]
                nom_args = {p.arg: p.nom for p in params if p not in ss}
            else:
                nom_args = {}
                ss = params

            lbmin, ubmax = float("inf"), -float("inf")
            for combo in product((min, max), repeat=len(ss)):
                kwargs = {p.arg: get(p.lb, p.ub) for get, p in zip(combo, ss)}
                result = func(*args, **kwargs, **f_kwargs, **nom_args)
                lbmin = result if result < lbmin else lbmin
                ubmax = result if result > ubmax else ubmax

            nom = func(*args, **{p.arg: p.nom for p in params}, **f_kwargs)
            lbmin = nom if nom < lbmin else lbmin
            ubmax = nom if nom > ubmax else ubmax

            return Param(arg=arg, nom=nom, lb=lbmin, ub=ubmax, **p_kwargs)

        return wrapper

    return decorator
