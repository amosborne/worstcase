from enum import Enum, auto
from inspect import Signature
from itertools import product

import networkx as nx  # type: ignore


class Mode(Enum):
    EV = auto()
    MC = auto()


class Param:
    def __init__(self, nom, lb, ub):
        self.nom = nom  # nominal value
        self.lb = lb  # lower bound
        self.ub = ub  # upper bound

    def build(self):
        return self


class ParamBuilder(Param):
    def __init__(self, func, mode, *args, **kwargs):
        sig = Signature.from_callable(func)

        self.func = func
        self.bind = sig.bind_partial(*args, **kwargs)
        self.mode = mode

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
    def byrange(nom, lb, ub):
        return Param(nom, lb, ub)

    @staticmethod
    def bytol(nom, tol, rel):
        tol = nom * tol if rel else tol
        return Param(nom, nom - tol, nom + tol)

    @classmethod
    def ev(cls, *args, **kwargs):
        def decorator(func):
            return cls(func, Mode.EV, *args, **kwargs)

        return decorator

    @classmethod
    def mc(cls, *args, **kwargs):
        def decorator(func):
            return cls(func, Mode.MC, *args, **kwargs)

        return decorator

    def __call__(self, *args, **kwargs):
        newsig = Signature.from_callable(self.func)
        newbind = newsig.bind_partial(*args, **kwargs)
        finalbind = {**self.bind.arguments}
        finalbind.update(newbind.arguments)

        try:
            bind = newsig.bind(**finalbind)
            params = list(
                filter(lambda v: isinstance(v, Param), bind.arguments.values())
            )
            if not params:
                return self.func(**finalbind)
        except TypeError:
            pass

        return ParamBuilder(self.func, self.mode, **finalbind)

    def build(self):
        predecessors = {
            k: v.build() for k, v in self.bind.arguments.items() if isinstance(v, Param)
        }

        if self.mode is Mode.EV:

            lbmin, ubmax = float("inf"), -float("inf")
            for combo in product((min, max), repeat=len(predecessors)):
                kwargs = {**self.bind.arguments}
                kwargs.update(
                    {
                        k: get(v.lb, v.ub)
                        for get, (k, v) in zip(combo, predecessors.items())
                    }
                )
                result = self.func(**kwargs)
                lbmin = result if result < lbmin else lbmin
                ubmax = result if result > ubmax else ubmax

            kwargs.update({k: v.nom for (k, v) in predecessors.items()})
            nom = self.func(**kwargs)
            return Param(nom=nom, lb=lbmin, ub=ubmax)

        else:
            pass
