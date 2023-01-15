from enum import Enum, auto
from inspect import Signature
from itertools import product

import networkx as nx
from pint import Quantity, UnitRegistry
from pyDOE import lhs

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
    def __init__(self, nom, lb, ub, tag, sigfig):
        nom = nom if isinstance(nom, Quantity) else nom * Unit([])
        lb = lb if isinstance(lb, Quantity) else lb * Unit([])
        ub = ub if isinstance(ub, Quantity) else ub * Unit([])

        assert lb.u == nom.u == ub.u, "Parameter bounds have inconsistent units."
        assert lb <= nom <= ub, "Parameter bounds are out of order."

        self.nom = nom  # nominal quantity
        self.lb = lb  # lower bound quantity
        self.ub = ub  # upper bound quantity
        self.tag = tag  # string identifier
        self.sigfig = sigfig  # string significant digits

    @staticmethod
    def byrange(nom, lb, ub, tag="", sigfig=4):
        return Parameter(nom, lb, ub, tag, sigfig)

    @staticmethod
    def bytol(nom, tol, rel, tag="", sigfig=4):
        tol = nom * tol if rel else tol
        return Parameter(nom, nom - tol, nom + tol, tag, sigfig)

    def __repr__(self):
        pretty = "0.{sigfig}G~P".format(sigfig=self.sigfig)
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

    def graph(self):
        graph = nx.DiGraph()
        graph.add_node(self, latest=self)
        return graph


class By(Enum):
    EV = auto()
    MC = auto()


class Derivative(AbstractParameter):
    def __init__(self, func, method, tag, sigfig, n, *args, **kwargs):
        self.func = func
        self.method = method
        self.tag = tag
        self.sigfig = sigfig
        self.n = n

        funcsig = Signature.from_callable(func)
        binding = funcsig.bind_partial(*args, **kwargs)
        self.kwargs = binding.arguments

    @staticmethod
    def byev(*args, tag="", sigfig=4, n=None, **kwargs):
        def decorator(func):
            return Derivative(func, By.EV, tag, sigfig, n, *args, **kwargs)

        return decorator

    @staticmethod
    def bymc(*args, tag="", sigfig=4, n=1000, **kwargs):
        def decorator(func):
            return Derivative(func, By.MC, tag, sigfig, n, *args, **kwargs)

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

    def __repr__(self):
        return self.derive().__repr__()

    def __call__(self, *args, tag=None, sigfig=None, n=None, **kwargs):

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
        sigfig = self.sigfig if sigfig is None else sigfig
        n = self.n if n is None else n
        return Derivative(self.func, self.method, tag, sigfig, n, **new_kwargs)

    def graph(self, ss=None):

        # Identify the AbstractParameter predecessors to this Derivative.
        is_abstractparam = lambda arg: isinstance(arg, AbstractParameter)  # noqa: E731
        preds = filter(is_abstractparam, self.kwargs.values())

        # Build the directed graph recursively, where:
        # - each node is an AbstractParameter
        # - each node has a "latest" field to hold the latest
        #   AbstractParameter assignment which will be used in evaluation
        # - each directed edge has a "method" field (By.MC or By.EV)
        graph = nx.DiGraph()
        graph.add_node(self, latest=self)
        for pred in preds:
            graph.add_edge(pred, self, method=self.method)
            graph = nx.compose(graph, pred.graph())

        # A sensitivity study is provided by the "ss" argument, a list of
        # AbstractParameters. A sensitivity study operates only over Parameters;
        # if Derivatives are provided they are used only to identify those
        # Parameters which are ancestors to those Derivatives. Any Parameters
        # which are not identified for the sensitivity study are replaced such
        # that their bounds are the same as their nominal value.
        if ss is not None:
            ss = set(ss) if isinstance(ss, list) else {ss}
            is_param = lambda arg: isinstance(arg, Parameter)  # noqa: E731
            ss_params = set()
            for abstractparam in ss:
                if is_param(abstractparam):
                    ss_params.add(abstractparam)
                else:
                    ancestors = nx.ancestors(graph, abstractparam)
                    ss_params |= set(filter(is_param, ancestors))

            for node in graph.nodes:
                if is_param(node) and node not in ss_params:
                    latest = Parameter(
                        node.nom, node.nom, node.nom, node.tag, node.sigfig
                    )
                    graph.nodes[node]["latest"] = latest

        return graph

    def ss(self, ss):
        return self.derive(ss)

    def derive(self, ss=None):

        # Construct a directed graph representing the interdependency of
        # all AbstractParameters needed to derive this Derivative.
        graph = self.graph(ss)
        cycles = list(nx.simple_cycles(graph))
        assert not cycles, "Derivative cannot have cyclical dependencies."

        # Traverse the graph (in any order). For each node, get the set
        # of ancestor nodes. If no ancestor node contains an out-edge
        # to a node not already in the set of ancestor nodes plus the
        # the current descendant node, flag the current descendant node.
        for node in graph.nodes:
            if isinstance(node, Parameter):
                graph.nodes[node]["flag"] = False
            else:
                ancestors = nx.ancestors(graph, node)
                out_edges = graph.out_edges(ancestors)
                out_nodes = {edge[1] for edge in out_edges}
                eval_group = ancestors.copy()
                eval_group.add(node)
                graph.nodes[node]["flag"] = eval_group >= out_nodes

        # Continuously loop through the nodes in the graph flagged for
        # evaluation until only the final derived Parameter remains.
        # Skip nodes depending on nodes flagged for evaluation.
        while not isinstance(graph.nodes[self]["latest"], Parameter):
            eval_nodes = [n for (n, d) in graph.nodes(data=True) if d["flag"]]
            for node in eval_nodes:
                ancestors = nx.ancestors(graph, node)
                if not any([anc in eval_nodes for anc in ancestors]):
                    # Evaluate using the current nodes method (EV or MC).
                    # This will overrule any depending node's method.
                    if node.method is By.EV:
                        param = extreme_value(graph, node)
                    else:
                        param = monte_carlo(graph, node)

                    # Update the current node with the latest evaluation.
                    graph.nodes[node]["latest"] = param
                    graph.nodes[node]["flag"] = False

        return graph.nodes[self]["latest"]


def eval_graph(graph, eval_node, eval_init):

    # Construct the evaluation group and initialize.
    graph.nodes[eval_node]["val"] = None
    ancestors = set(nx.ancestors(graph, eval_node)) | {eval_node}
    eval_group = set()
    for anc in ancestors:
        if isinstance(graph.nodes[anc]["latest"], Derivative):
            graph.nodes[anc]["val"] = None
            eval_group.add(anc)
            for pred in graph.predecessors(anc):
                if isinstance(graph.nodes[pred]["latest"], Parameter):
                    graph.nodes[pred]["val"] = eval_init[pred]
                    eval_group.add(pred)

    # Continuously bubble-up evaluations to the eval-node.
    while graph.nodes[eval_node]["val"] is None:
        for node in eval_group:
            # Skip any nodes already evaluated.
            if not graph.nodes[node]["val"] is None:
                continue

            # Skip any nodes with unevaluated predecessors.
            if any(
                graph.nodes[pred]["val"] is None for pred in graph.predecessors(node)
            ):
                continue

            # Construct the keyword arguments.
            kwargs = node.kwargs.copy()
            for k, v in kwargs.items():
                if v in graph.predecessors(node):
                    kwargs[k] = graph.nodes[v]["val"]

            graph.nodes[node]["val"] = node.func(**kwargs)

    # Return the eval-node final evaluation.
    return graph.nodes[eval_node]["val"]


def eval_nominal(graph, eval_node):
    # Get the Parameters to initialize of the eval-node.
    ancestors = set(nx.ancestors(graph, eval_node)) | {eval_node}
    params = set()
    for anc in ancestors:
        if isinstance(graph.nodes[anc]["latest"], Derivative):
            for pred in graph.predecessors(anc):
                if isinstance(graph.nodes[pred]["latest"], Parameter):
                    params.add(pred)

    # Calculate the nominal value.
    eval_init = {p: p.nom for p in params}
    nom = eval_graph(graph, eval_node, eval_init)

    return params, nom


def extreme_value(graph, eval_node):
    params, nom = eval_nominal(graph, eval_node)

    # Loop through all max/min combinations for all primitives.
    lbmin, ubmax = float("inf") * nom.u, -float("inf") * nom.u
    for combo in product((min, max), repeat=len(params)):
        eval_init = {}
        for p, c in zip(params, combo):
            latest = graph.nodes[p]["latest"]
            eval_init[p] = c(latest.lb, latest.ub)
        result = eval_graph(graph, eval_node, eval_init)
        lbmin = result if result < lbmin else lbmin
        ubmax = result if result > ubmax else ubmax

    return Parameter(nom, lbmin, ubmax, eval_node.tag, eval_node.sigfig)


def monte_carlo(graph, eval_node):
    params, nom = eval_nominal(graph, eval_node)

    # Run n Monte Carlo evaluations with Latin Hypercube Sampling.
    matrix = lhs(len(params), samples=eval_node.n)
    results = []
    for row in matrix:
        eval_init = {}
        for p, x in zip(params, row):
            latest = graph.nodes[p]["latest"]
            eval_init[p] = x * (latest.ub - latest.lb) + latest.lb
        results.append(eval_graph(graph, eval_node, eval_init))

    return Parameter(nom, min(results), max(results), eval_node.tag, eval_node.sigfig)
