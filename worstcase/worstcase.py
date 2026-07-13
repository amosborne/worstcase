from collections import namedtuple
from enum import Enum, auto
from inspect import Signature
from itertools import product
from warnings import warn

import networkx as nx
import numpy as np
from pydoe import lhs

from worstcase.unit_adapter import unit_adapter


class AbstractParameter:
    @property
    def units(self):
        """Return the unit of the Parameter"""
        return unit_adapter(type(self.nom)).units(self.nom)

    @property
    def u(self):
        """Return the unit of the Parameter"""
        return self.units

    @staticmethod
    def get_units(value):
        """Return the unit of a value."""
        return unit_adapter(type(value)).units(value)

    @staticmethod
    def get_val(value):
        """Return the unitless value of a dimensioned value"""
        return unit_adapter(type(value)).magnitude(value)

    @staticmethod
    def compare_units(val1, val2):
        """Check if two values have compatible units.

        Arguments:
            val1, val2: values to compare

        Returns:
            True if the values have compatible units
            False if the values have incompatible units
        """
        return unit_adapter(type(val1)).compatible(val1, val2)

    @staticmethod
    def convert_units(val1, val2):
        """Convert val1 to the units of val2"""
        return unit_adapter(type(val1)).convert(val1, val2)

    @staticmethod
    def is_dimensionless(val):
        """Check if a value is dimensionless"""
        return unit_adapter(type(val)).is_dimensionless(val)

    def equivalent(self, other):
        """Check if one Parameter is equivalent to another.

        Equivalency means that nom, lb, and ub are the same,
            while the tag, derivation, and sigfigs may differ.

        Returns:
            boolean; True if they are equivalent, false otherwise
        """
        if not isinstance(other, AbstractParameter):
            return False
        if type(self.units) is not type(other.units):
            return False
        return (
            (self.nom == other.nom) and (self.lb == other.lb) and (self.ub == other.ub)
        )


Derivation = namedtuple("Derivation", "nom lb ub")


class Parameter(AbstractParameter):
    def __init__(
        self, nom, lb, ub, tag, sigfig, derivation=Derivation(nom={}, lb={}, ub={})
    ):
        """Create a new Parameter.

        Parameters:
            nom: nominal value of the parameter
            lb: lower bound of the parameter
            ub: upper bound of the parameter
            tag: display identifier for the parameter
            sigfig: number of significant figures to display as a string
            derivation: values used in derivation of this value

        Raises:
            ValueError: if the nom/lb/ub values are inconsistent
        """
        if not (
            AbstractParameter.compare_units(lb, nom)
            and AbstractParameter.compare_units(ub, nom)
        ):
            raise ValueError("Parameter bounds have inconsistent units.")
        if not (lb <= nom <= ub):
            raise ValueError("Parameter bounds are out of order.")

        self.nom = nom
        self.lb = lb
        self.ub = ub
        self.tag = tag
        self.sigfig = sigfig
        self.derivation = derivation

    @staticmethod
    def byrange(nom, lb, ub, tag="", sigfig=4):
        """Define a parameter by an upper and lower bound.

        Parameters:
            nom: nominal value
            lb: lower bound on the nominal value
            ub: upper bound on the nominal value
            tag: name of the parameter
            sigfig: number of significant figures to print
        """
        if not (
            AbstractParameter.compare_units(nom, lb)
            and AbstractParameter.compare_units(nom, ub)
        ):
            raise ValueError("Nominal value and bounds must have compatible units")
        lb = AbstractParameter.convert_units(lb, nom)
        ub = AbstractParameter.convert_units(ub, nom)
        return Parameter(nom, lb, ub, tag, sigfig)

    @staticmethod
    def bytol(nom, tol, rel, tag="", sigfig=4):
        """Define a parameter by a tolerance.

        Aliases `bytol_rel` and `bytol_abs` with boolean to pick between the two.

        Parameters:
            nom: nominal value
            tol: tolerance to be applied to nominal value
            rel: boolean, defines if tol is a relative value (rel=true)
                or an absolute value (rel=false)
            tag: name of the parameter
            sigfig: number of significant figures to print
        """
        if rel:
            return Parameter.bytol_rel(nom, tol, tag, sigfig)
        else:
            return Parameter.bytol_abs(nom, tol, tag, sigfig)

    @staticmethod
    def bytol_rel(nom, rel_tol, tag="", sigfig=4):
        """Define a parameter by a relative tolerance.

        ex. bytol_rel(5, 0.1) would result in a nom=5, lb=4.5, ub=5.5.

        Parameters:
            nom: nominal value
            rel_tol: relative tolerance to be applied to nominal value
            tag: name of the parameter
            sigfig: number of significant figures to print
        """
        if not AbstractParameter.is_dimensionless(rel_tol):
            raise ValueError("Relative value cannot have units")

        tol = nom * rel_tol
        return Parameter(nom, nom - tol, nom + tol, tag, sigfig)

    @staticmethod
    def bytol_abs(nom, abs_tol, tag="", sigfig=4):
        """Define a parameter by an absolute tolerance.

        ex. bytol_abs(5, 0.1) would result in a nom=5, lb=4.9, ub=5.1.

        Parameters:
            nom: nominal value
            abs_tol: absolute tolerance to be applied to nominal value
            tag: name of the parameter
            sigfig: number of significant figures to print
        """

        if not AbstractParameter.compare_units(nom, abs_tol):
            raise ValueError("Nominal and absolute values must have compatible units")
        abs_tol = AbstractParameter.convert_units(abs_tol, nom)
        return Parameter(nom, nom - abs_tol, nom + abs_tol, tag, sigfig)

    def formatter_string(self):
        """Generate the formatter string for the current data"""
        return unit_adapter(type(self.nom)).format_str(self.sigfig)

    def __repr__(self):
        pretty = self.formatter_string()
        tag = "{tag}: ".format(tag=self.tag) if self.tag else ""
        nom = ("{nom:" + pretty + "} (nom), ").format(nom=self.nom)
        lb = ("{lb:" + pretty + "} (lb), ").format(lb=self.lb)
        ub = ("{ub:" + pretty + "} (ub)").format(ub=self.ub)

        return tag + nom + lb + ub

    def _repr_latex_(self):
        upper_deviation = self.ub - self.nom
        lower_deviation = self.nom - self.lb
        pretty = self.formatter_string()
        if f"{upper_deviation:{pretty}}" == f"{lower_deviation:{pretty}}":
            return f"${self.nom:{pretty}} \\pm {upper_deviation:{pretty}}$"
        else:
            return (
                f"${self.nom:{pretty}} \\,{{}}"
                f"^{{+{upper_deviation:{pretty}}}}"
                f"_{{-{lower_deviation:{pretty}}}}$"
            )

    def __format__(self, fmt_spec):
        if not fmt_spec:
            return str(self)
        if fmt_spec == "L":
            return self._repr_latex_()
        else:
            raise ValueError("Invalid format specifier")

    def apply(self, funcname, *args, **kwargs):
        """Apply a function to nominal and upper/lower bound values.

        Ex. `self.apply(ito_base_units)` will convert all values to base units
            if using Pint values.
        """
        getattr(self.nom, funcname)(*args, **kwargs)
        getattr(self.lb, funcname)(*args, **kwargs)
        getattr(self.ub, funcname)(*args, **kwargs)
        return self

    def graph(self):
        graph = nx.DiGraph()
        graph.add_node(self, latest=self)
        return graph


class By(Enum):
    EV = auto()
    MC = auto()
    RSS = auto()


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

    @staticmethod
    def byrss(*args, tag="", sigfig=4, n=None, **kwargs):
        def decorator(func):
            return Derivative(func, By.RSS, tag, sigfig, n, *args, **kwargs)

        return decorator

    @property
    def nom(self):
        return eval_nominal(self.graph(), self)[1]

    @property
    def lb(self):
        return self.derive().lb

    @property
    def ub(self):
        return self.derive().ub

    @property
    def derivation(self):
        return self.derive().derivation

    def __repr__(self):
        return self.derive().__repr__()

    def _repr_latex_(self):
        return self.derive()._repr_latex_()

    def __format__(self, format_spec):
        return self.derive().__format__(format_spec)

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
        graph = nx.DiGraph()
        graph.add_node(self, latest=self)
        for pred in preds:
            graph.add_edge(pred, self)
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
        if cycles:
            raise RecursionError("Derivative cannot have cyclical dependencies.")

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
                    # Evaluate using the current nodes method (EV, MC, or RSS).
                    # This will overrule any depending node's method.
                    if node.method is By.EV:
                        param = extreme_value(graph, node)
                    elif node.method is By.MC:
                        param = monte_carlo(graph, node)
                    else:
                        param = root_sum_square(graph, node)

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

    derivation_nom = {p: p.nom for p in params}

    # Loop through all max/min combinations for all primitives.
    nom_units = AbstractParameter.get_units(nom)
    lbmin, ubmax = float("inf") * nom_units, -float("inf") * nom_units
    for combo in product((min, max), repeat=len(params)):
        eval_init = {}
        for p, c in zip(params, combo):
            latest = graph.nodes[p]["latest"]
            eval_init[p] = c(latest.lb, latest.ub)
        result = eval_graph(graph, eval_node, eval_init)

        if result < lbmin:
            lbmin = result
            derivation_lb = eval_init

        if result > ubmax:
            ubmax = result
            derivation_ub = eval_init

    derivation = Derivation(nom=derivation_nom, lb=derivation_lb, ub=derivation_ub)

    return Parameter(nom, lbmin, ubmax, eval_node.tag, eval_node.sigfig, derivation)


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

    # todo: include derivation in the resulting parameter
    return Parameter(nom, min(results), max(results), eval_node.tag, eval_node.sigfig)


def root_sum_square(graph, eval_node):
    params, nom = eval_nominal(graph, eval_node)

    # For each parameter, calculate the result using both the lower and upper
    # bounds while holding all other parameters at their nominal value.
    lbs = []
    ubs = []
    for pvaried in params:
        eval_init = {p: graph.nodes[p]["latest"].nom for p in params - {pvaried}}
        eval_init[pvaried] = graph.nodes[pvaried]["latest"].lb
        ret_lb = AbstractParameter.get_val(eval_graph(graph, eval_node, eval_init))
        eval_init[pvaried] = graph.nodes[pvaried]["latest"].ub
        ret_ub = AbstractParameter.get_val(eval_graph(graph, eval_node, eval_init))
        ubs.append(max([ret_lb, ret_ub]))
        lbs.append(min([ret_lb, ret_ub]))

    # Define a metric to detect an asymmetric result and warn the user.
    nom_val = AbstractParameter.get_val(nom)
    nom_units = AbstractParameter.get_units(nom)
    ub_deviation = np.array(ubs) - nom_val
    lb_deviation = nom_val - np.array(lbs)
    mean_deviation = (ub_deviation + lb_deviation) / 2
    max_deviation = np.maximum(ub_deviation, lb_deviation)
    if np.any(max_deviation > mean_deviation * 1.05):
        warn("RSS method is not recommended for asymmetric derived parameters.")

    # Compute the RSS (root-sum-square).
    ub = nom_val + np.sqrt(np.sum(np.square(np.array(ubs) - nom_val)))
    lb = nom_val - np.sqrt(np.sum(np.square(np.array(lbs) - nom_val)))

    # todo: include derivation in the resulting parameter
    return Parameter(
        nom,
        lb * nom_units,
        ub * nom_units,
        eval_node.tag,
        eval_node.sigfig,
    )
