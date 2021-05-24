from enum import Enum, auto
from inspect import Signature
from itertools import product
from types import SimpleNamespace

import networkx as nx  # type: ignore
from pint import Quantity, UnitRegistry  # type: ignore
from pyDOE import lhs  # type: ignore
from treelib import Tree  # type: ignore

Config = SimpleNamespace(n=5000, sigfig=3)
Unit = UnitRegistry()


class Mode(Enum):
    EV = auto()
    MC = auto()


class Param:
    def __init__(self, nom, lb, ub, tag):
        self.nom = nom  # nominal value
        self.lb = lb  # lower bound
        self.ub = ub  # upper bound
        self.tag = tag  # string identifier

    def build(self):
        return self

    def __repr__(self):
        pretty = "0.{sigfig}g~P".format(sigfig=Config.sigfig)
        tag = "{tag}: ".format(tag=self.tag) if self.tag else ""
        nom = ("{nom:" + pretty + "} (nom), ").format(nom=self.nom.to_compact())
        lb = ("{lb:" + pretty + "} (min), ").format(lb=self.lb.to_compact())
        ub = ("{ub:" + pretty + "} (max)").format(ub=self.ub.to_compact())
        return tag + nom + lb + ub


class ParamBuilder(Param):
    def __init__(self, func, mode, *args, tag="", **kwargs):
        sig = Signature.from_callable(func)

        self.func = func
        self.bind = sig.bind_partial(*args, **kwargs)
        self.mode = mode
        self.tag = tag

        predecessors = filter(
            lambda v: isinstance(v, Param), self.bind.arguments.values()
        )

        graph = nx.DiGraph()
        graph.add_node(self)
        for pred in predecessors:
            graph.add_node(pred)
            graph.add_edge(pred, self, mode=mode)
            if isinstance(pred, ParamBuilder):
                graph = nx.compose(graph, pred.graph)

        assert nx.is_tree(graph), "Composed worst case analysis must be acyclic."
        self.graph = graph

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
            newtag = func.__name__ if tag is None else tag
            return cls(func, Mode.EV, *args, **kwargs, tag=newtag)

        return decorator

    @classmethod
    def mc(cls, *args, tag=None, **kwargs):
        def decorator(func):
            newtag = func.__name__ if tag is None else tag
            return cls(func, Mode.MC, *args, **kwargs, tag=newtag)

        return decorator

    def __call__(self, *args, tag=None, **kwargs):
        # If no arguments, return the built parameter.
        if not args and not kwargs:
            return self.build()

        # Otherwise, update the binding arguments.
        newsig = Signature.from_callable(self.func)
        newbind = newsig.bind_partial(*args, **kwargs)
        finalbind = {**self.bind.arguments}
        finalbind.update(newbind.arguments)

        # If all arguments are not parameters, return a single value.
        try:
            bind = newsig.bind(**finalbind)
            params = [p for p in bind.arguments.values() if isinstance(p, Param)]
            if not params:
                return self.func(**finalbind)
        except TypeError:
            pass

        # Otherwise, return a new parameter builder.
        tag = self.func.__name__ if tag is None else tag
        return ParamBuilder(self.func, self.mode, tag=tag, **finalbind)

    def build(self):
        predecessors = {
            k: v.build() for k, v in self.bind.arguments.items() if isinstance(v, Param)
        }

        kwargs = {**self.bind.arguments}
        kwargs.update({k: v.nom for (k, v) in predecessors.items()})
        nom = self.func(**kwargs)

        # EXTREME VALUE ANALYSIS
        if self.mode is Mode.EV:
            lbmin, ubmax = float("inf"), -float("inf")
            for combo in product((min, max), repeat=len(predecessors)):
                kwargs.update(
                    {
                        k: get(v.lb, v.ub)
                        for get, (k, v) in zip(combo, predecessors.items())
                    }
                )
                result = self.func(**kwargs)
                if not isinstance(lbmin, Quantity):
                    lbmin *= result.units
                    ubmax *= result.units

                lbmin = result if result < lbmin else lbmin
                ubmax = result if result > ubmax else ubmax

            return Param(nom=nom, lb=lbmin, ub=ubmax, tag=self.tag)

        # MONTE CARLO ANALYSIS
        else:
            matrix = lhs(len(predecessors), samples=Config.n)
            results = []

            for row in matrix:
                kwargs.update(
                    {
                        k: (x * (v.ub - v.lb) + v.lb)
                        for (x, (k, v)) in zip(row, predecessors.items())
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
