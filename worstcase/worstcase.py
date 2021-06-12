from enum import Enum, auto
from inspect import Signature
from itertools import product
from types import SimpleNamespace

import networkx as nx  # type: ignore
from pint import Quantity, UnitRegistry  # type: ignore
from pyDOE import lhs  # type: ignore
from treelib import Tree  # type: ignore

Config = SimpleNamespace(n=5000, sigfig=4)
Unit = UnitRegistry()


class AbstractParameter:
    def check(self, dimensionality):
        return self.nom.check(dimensionality)

    @property
    def dimensionality(self):
        return self.nom.dimensionality

    def is_compatible_with(self, other, *contexts, **ctx_kwargs):
        return self.nom.is_compatible_with(other, *contexts, **ctx_kwargs)

    @property
    def units(self):
        return self.nom.units

    @property
    def u(self):
        return self.units


class Parameter(AbstractParameter):
    def __init__(self, nom, lb, ub, tag):
        nom = nom if isinstance(nom, Quantity) else nom * Unit([])
        lb = lb if isinstance(lb, Quantity) else lb * Unit([])
        ub = ub if isinstance(ub, Quantity) else ub * Unit([])

        assert lb.u == nom.u == ub.u, "Parameter bounds have inconsistent units."
        assert lb <= nom <= ub, "Parameter bounds are out of order."

        self.nom = nom  # nominal quantity
        self.lb = lb  # lower bound quantity
        self.ub = ub  # upper bound quantity
        self.tag = tag  # string identifier

    @staticmethod
    def byrange(nom, lb, ub, tag=""):
        return Parameter(nom, lb, ub, tag)

    @staticmethod
    def bytol(nom, tol, rel, tag=""):
        tol = nom * tol if rel else tol
        return Parameter(nom, nom - tol, nom + tol, tag)

    def __repr__(self):
        pretty = "0.{sigfig}G~P".format(sigfig=Config.sigfig)
        tag = "{tag}: ".format(tag=self.tag) if self.tag else ""
        nom = ("{nom:" + pretty + "} (nom), ").format(nom=self.nom.to_compact())
        lb = ("{lb:" + pretty + "} (lb), ").format(lb=self.lb.to_compact())
        ub = ("{ub:" + pretty + "} (ub)").format(ub=self.ub.to_compact())
        return tag + nom + lb + ub

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


class By(Enum):
    EV = auto()
    MC = auto()


class Derivative(AbstractParameter):
    def __init__(self, func, method, tag, *args, **kwargs):
        self.func = func
        self.method = method
        self.tag = tag

        funcsig = Signature.from_callable(func)
        binding = funcsig.bind_partial(*args, **kwargs)
        self.kwargs = binding.arguments

    @staticmethod
    def byev(*args, tag="", **kwargs):
        def decorator(func):
            return Derivative(func, By.EV, tag, *args, **kwargs)

        return decorator

    @staticmethod
    def bymc(*args, tag="", **kwargs):
        def decorator(func):
            return Derivative(func, By.MC, tag, *args, **kwargs)

        return decorator

    @property
    def nom(self):
        return self.derive().nom

    @property
    def lb(self):
        return self.derive().lb

    @property
    def ub(self):
        return self.derive().ub

    def __call__(self, *args, tag=None, **kwargs):
        # If args/kwargs form a complete binding with no AbstractParameters
        # then return a simple call to the underlying function.
        try:
            funcsig = Signature.from_callable(self.func)
            binding = funcsig.bind(*args, **kwargs)
            kwargs = binding.arguments

            is_abstract_parameter = [
                isinstance(v, AbstractParameter) for v in kwargs.values()
            ]
            if not any(is_abstract_parameter):
                return self.func(**kwargs)
        except TypeError:
            pass

        # If no arguments are supplied then return the derived Parameter.
        if not any([args, kwargs, tag]):
            return self.derive()

        # Otherwise return a new Derivative with an updated binding.
        new_kwargs = self.kwargs.copy()
        funcsig = Signature.from_callable(self.func)
        binding = funcsig.bind_partial(*args, **kwargs)
        new_kwargs.update(binding.arguments)

        tag = self.tag if tag is None else tag
        return Derivative(self.func, self.method, tag, **new_kwargs)

    def derive(self):
        # Construct a directed graph representing the interdependency of
        # all AbstractParameters needed to derive this Derivative.
        graph = self.graph
        cycles = list(nx.simple_cycles(graph))
        assert not cycles, "Derivative cannot have cyclical dependencies."
        nx.draw(graph)

        # Traverse the graph (in any order). For each node, get the set
        # of ancestor nodes. If no ancestor node contains an out-edge
        # to a node not already in the set of ancestor nodes plus the
        # the current descendant node, flag the current descendant node.
        for node in graph.nodes:
            if isinstance(node, Parameter):
                graph.nodes[node]["eval"] = False
                continue

            ancestors = nx.ancestors(graph, node)
            eval_group = ancestors.copy()
            eval_group.add(node)
            out_edges = graph.out_edges(ancestors)
            out_nodes = {e[1] for e in out_edges}

            graph.nodes[node]["eval"] = eval_group >= out_nodes

        # Continuously loop through the nodes in the graph flagged for
        # evaluation until only the final derived Parameter remains.
        while len(graph.nodes) > 1:
            eval_nodes = [n for (n, d) in graph.nodes(data=True) if d["eval"]]

            for node in eval_nodes:
                # Skip nodes depending on nodes flagged for evaluation.
                ancestors = nx.ancestors(graph, node)
                if any([anc in eval_nodes for anc in ancestors]):
                    continue

                # Evaluate using the current nodes method (EV or MC).
                # This will overrule any depending node's method.
                if node.method is By.EV:
                    param = extreme_value(graph, node)
                else:
                    param = monte_carlo(graph, node)

                # Replace the ancestor nodes and current node with
                # the newly derived Parameter.
                graph.remove_nodes_from(ancestors)
                graph.remove_node(node)
                graph.add_node(param, **{"eval": False})

        return list(graph.nodes)[0]

    @property
    def graph(self):
        graph = nx.DiGraph()
        graph.add_node(self)
        predecessors = [
            v for v in self.kwargs.values() if isinstance(v, AbstractParameter)
        ]
        for predecessor in predecessors:
            graph.add_node(predecessor)
            graph.add_edge(predecessor, self, method=self.method)
            if isinstance(predecessor, Derivative):
                graph = nx.compose(graph, predecessor.graph)

        return graph


def extreme_value(graph, node):
    return 1


def monte_carlo(graph, node):
    return None
