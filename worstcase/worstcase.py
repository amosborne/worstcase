from enum import Enum, auto
from inspect import Signature

import networkx as nx  # type: ignore


class Mode(Enum):
    EV = auto()
    MC = auto()


class Param:
    def __init__(self, nom, lb, ub):
        self.nom = nom  # nominal value
        self.lb = lb  # lower bound
        self.ub = ub  # upper bound

        graph = nx.DiGraph()
        graph.add_node(self)
        self.graph = graph


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
            graph = nx.compose(graph, pred.graph)

        # VERIFY: No node is used in both MC and EV worst-case analyses.
        for node in graph:
            out_edges = graph.out_edges(node, data="mode")
            if out_edges:
                modes = [m for (_, _, m) in out_edges]
                assert modes.count(modes[0]) == len(
                    modes
                ), "A parameter cannot be used in both MC and EV analysis modes."

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
        return ParamBuilder(self.func, self.mode, **finalbind)

    def build(self):
        raise NotImplementedError("Param building not implemented.")
