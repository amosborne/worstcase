from enum import Enum, auto
from inspect import Signature
from itertools import product
from types import SimpleNamespace

import networkx as nx  # type: ignore
from pint import UnitRegistry  # type: ignore
from pyDOE import lhs  # type: ignore
from treelib import Tree  # type: ignore

Config = SimpleNamespace(n=5000, sigfig=3)
Unit = UnitRegistry()


class Mode(Enum):
    EV = auto()
    MC = auto()


class Param:
    def __init__(self, nom, lb, ub, tag):
        assert lb <= nom <= ub, "Parameter bounds out of order."
        assert lb.u == nom.u == ub.u

        self.nom = nom  # nominal value
        self.lb = lb  # lower bound
        self.ub = ub  # upper bound
        self.tag = tag  # string identifier

    def __repr__(self):
        pretty = "0.{sigfig}g~P".format(sigfig=Config.sigfig)
        tag = "{tag}: ".format(tag=self.tag) if self.tag else ""
        nom = ("{nom:" + pretty + "} (nom), ").format(nom=self.nom.to_compact())
        lb = ("{lb:" + pretty + "} (lb), ").format(lb=self.lb.to_compact())
        ub = ("{ub:" + pretty + "} (ub)").format(ub=self.ub.to_compact())
        return tag + nom + lb + ub

    def check(self, dimension):
        return self.nom.check(dimension)

    @property
    def dimensionality(self):
        return self.nom.dimensionality

    def is_compatible_with(self, other, *contexts, **ctx_kwargs):
        return self.nom.is_compatible_with(other, *contexts, **ctx_kwargs)

    def ito(self, other=None, *contexts, **ctx_kwargs):
        self.lb.ito(other, *contexts, **ctx_kwargs)
        self.nom.ito(other, *contexts, **ctx_kwargs)
        self.ub.ito(other, *contexts, **ctx_kwargs)
        return self

    def ito_base_units(self):
        self.lb.ito_base_units()
        self.nom.ito_base_units()
        self.ub.ito_base_units()
        return self

    def ito_reduced_units(self):
        self.lb.ito_reduced_units()
        self.nom.ito_reduced_units()
        self.ub.ito_reduced_units()
        return self

    def ito_root_units(self):
        self.lb.ito_root_units()
        self.nom.ito_root_units()
        self.ub.ito_root_units()
        return self

    @property
    def units(self):
        return self.nom.units

    @property
    def u(self):
        return self.units


def bind_kwargs(func, *args, **kwargs):
    sig = Signature.from_callable(func)
    binding = sig.bind_partial(*args, **kwargs)
    binding.apply_defaults()
    return binding.arguments


class ParamBuilder(Param):
    def __init__(self, func, mode, *args, tag=None, **kwargs):
        self.func = func
        self.mode = mode
        self.tag = func.__name__ if tag is None else tag
        self.kwargs = bind_kwargs(func, *args, **kwargs)

        assert nx.is_tree(self.graph), "Composed worst case analysis must be acyclic."

    def __call__(self, *args, tag=None, **kwargs):
        new_kwargs = self.kwargs.copy()
        new_kwargs.update(bind_kwargs(self.func, *args, **kwargs))

        # If a complete binding is made and none are params, return func eval
        if not [v for v in new_kwargs.values() if isinstance(v, Param)]:
            try:
                Signature.from_callable(self.func).bind(**new_kwargs)
                return self.func(**new_kwargs)
            except TypeError:
                pass

        # If no arguments are supplied, return the built param
        if not args and not kwargs and tag is None:
            return self.build()

        # Otherwise, return a new parambuilder with an updated binding
        tag = self.tag if tag is None else tag
        return ParamBuilder(self.func, self.mode, tag=tag, **new_kwargs)

    def ss(self, params, tag=None):
        params = [params] if not isinstance(params, list) else params

        kwargs = self.kwargs.copy()
        for k, v in kwargs.items():
            if v in params or not isinstance(v, Param):
                continue
            elif not isinstance(v, ParamBuilder):
                kwargs[k] = v.nom
            else:
                kwargs[k] = v.ss(params)

        tag = self.tag if tag is None else tag
        return self(tag=tag, **kwargs)

    @property
    def graph(self):
        graph = nx.DiGraph()
        graph.add_node(self)
        for pred in [v for v in self.kwargs.values() if isinstance(v, Param)]:
            graph.add_node(pred)
            graph.add_edge(pred, self, mode=self.mode)
            if isinstance(pred, ParamBuilder):
                graph = nx.compose(graph, pred.graph)

        return graph

    @property
    def nom(self):
        return self.build().nom

    @property
    def lb(self):
        return self.build().lb

    @property
    def ub(self):
        return self.build().ub

    @staticmethod
    def byrange(nom, lb, ub, tag="", unit=Unit([])):
        return Param(nom * unit, lb * unit, ub * unit, tag)

    @staticmethod
    def bytol(nom, tol, rel, tag="", unit=Unit([])):
        tol = nom * tol if rel else tol
        return Param(nom * unit, (nom - tol) * unit, (nom + tol) * unit, tag)

    @classmethod
    def ev(cls, *args, tag=None, **kwargs):
        def decorator(func):
            return cls(func, Mode.EV, *args, **kwargs, tag=tag)

        return decorator

    @classmethod
    def mc(cls, *args, tag=None, **kwargs):
        def decorator(func):
            return cls(func, Mode.MC, *args, **kwargs, tag=tag)

        return decorator

    def build(self):
        preds = {k: v for k, v in self.kwargs.items() if isinstance(v, Param)}
        preds_lbub = {k: (v.lb, v.ub) for k, v in preds.items()}
        preds_nom = {k: v.nom for k, v in preds.items()}

        kwargs = self.kwargs.copy()
        kwargs.update(preds_nom)
        nom = self.func(**kwargs)

        # EXTREME VALUE ANALYSIS
        if self.mode is Mode.EV:
            lbmin, ubmax = float("inf") * nom.u, -float("inf") * nom.u
            for combo in product((min, max), repeat=len(preds)):
                kwargs.update(
                    {k: get(*lbub) for get, (k, lbub) in zip(combo, preds_lbub.items())}
                )
                result = self.func(**kwargs)

                lbmin = result if result.m < lbmin.m else lbmin
                ubmax = result if result.m > ubmax.m else ubmax

            return Param(nom=nom, lb=lbmin, ub=ubmax, tag=self.tag)

        # MONTE CARLO ANALYSIS
        else:
            matrix = lhs(len(preds), samples=Config.n)
            results = []

            for row in matrix:
                kwargs.update(
                    {
                        k: (x * (ub - lb) + lb)
                        for (x, (k, (lb, ub))) in zip(row, preds_lbub.items())
                    }
                )
                results.append(self.func(**kwargs))

            return Param(nom=nom, lb=min(results), ub=max(results), tag=self.tag)

    def _make_tree(self):
        tree = Tree()
        mode = " (mc)" if self.mode is Mode.MC else " (ev)"
        tree.create_node(self.tag + mode, hash(self))

        for (u, v) in self.graph.in_edges(self):
            if isinstance(u, ParamBuilder):
                tree.paste(hash(self), u._make_tree())
            else:
                tree.create_node(u.tag, hash(u), parent=hash(self))

        return tree

    def __repr__(self):
        return str(self._make_tree())
