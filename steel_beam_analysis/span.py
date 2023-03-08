# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
from operator import attrgetter
from sortedcontainers import SortedSet

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis import units
import steel_beam_analysis.stringFixer as sf
from steel_beam_analysis.node import Node
from steel_beam_analysis.element import Element


class Span:
    """Model a simple span where the i node is the left support and the j node
    is the right support. Or, model a cantilever span where the i node or j node
    is a free node, and the other node is a support."""

    def __init__(self, iNode, jNode):
        if (type(iNode) != Node) or (type(jNode) != Node):
            raise TypeError

        self.cantilever = None
        self.deflChecks = {}
        self.deflLimits = None
        self.deflRatioLimits = {}
        self.deflRatios = {}
        self._elements = SortedSet([])
        self.iNode = iNode
        self.jNode = jNode
        self.L = self.jNode.loc - self.iNode.loc
        self.maxDeflNodes = {}
        self._nodes = None
        self.nodes = SortedSet([])
        self.ultDefls = None
        self.checkCantilever()

    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, val):
        if not isinstance(val, Element):
            raise TypeError
        self._elements = val

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, vals):
        if not isinstance(vals, SortedSet) and not isinstance(vals, list):
            raise TypeError
        if any(not isinstance(node, Node) for node in vals):
            raise TypeError
        self._nodes = vals

    def checkDeflections(self, type):
        """Check that deflections in the span are within limits."""

        ratio = self.deflRatios[type].to_reduced_units().magnitude
        if ratio:
            if ratio >= self.deflRatioLimits[type]:
                self.comparison = "<"
                result = "OK"
            else:
                self.comparison = (
                    ">"  # sign flipped because comparison done with reciprocal
                )
                result = "NG"
        else:
            self.comparison = "N/A"
            result = "Not calculated"
            result = "OK"

        self.deflChecks[type] = result

        if type == "LL":
            maxDefl = round(self.maxDeflNodes["LL"].deflMaxAbsL["val"], 2)
        elif type == "TL":
            maxDefl = round(self.maxDeflNodes["TL"].deflMaxAbs["val"], 2)
        equation = f"\\Delta_{{max}} = {maxDefl} = \\cfrac{{L}}{{{int(ratio)}}} {self.comparison} \\cfrac{{L}}{{{self.deflRatioLimits[type]}}} _sp_ _bf_{{({result})}}"
        return f"{type.title()} Deflection Check: \r${sf.fixUnits(equation)}$"

    def isDeflectionOK(self):
        return "NG" if "NG" in self.deflChecks.values() else "OK"

    def setTLdeflRatios(self, ratioLimit):
        """Set the total load deflection ratio and limiting ratio in the span."""
        self.deflRatios["TL"] = abs(self.L / self.maxDeflNodes["TL"].deflMaxAbs["val"])
        self.deflRatioLimits["TL"] = ratioLimit
        self.deflRatioLimits["TL"] /= (
            2 if self.cantilever else self.deflRatioLimits["TL"]
        )

    def setLLdeflRatios(self, ratioLimit):
        """Set the live load deflection ratio and limiting ratio in the span."""
        self.deflRatios["LL"] = abs(self.L / self.maxDeflNodes["LL"].deflMaxAbsL["val"])
        self.deflRatioLimits["LL"] = ratioLimit
        self.deflRatioLimits["LL"] /= (
            2 if self.cantilever else self.deflRatioLimits["LL"]
        )

    def setMaxTLdeflNode(self):
        """Set the node with the max total load deflection in the span."""
        self.maxDeflNodes["TL"] = max(
            self.nodes, key=lambda node: abs(node.deflMaxAbs["val"])
        )

    def setMaxLLdeflNode(self):
        """Set the node with the max live load deflection in the span."""
        self.maxDeflNodes["LL"] = max(
            self.nodes, key=lambda node: abs(node.deflMaxAbsL["val"])
        )

    def checkCantilever(self):
        """If either the i node or the j node of the span is free, report that
        the span is a cantilever."""
        if self.iNode.condition == "free" or self.jNode.condition == "free":
            self.cantilever = True
        else:
            self.cantilever = False

    def __cmp__(self, other):
        return cmp(self.iNode.loc, other.iNode.loc)

    def __lt__(self, other):
        return self.iNode.loc < other.iNode.loc

    def __gt__(self, other):
        return self.iNode.loc > other.iNode.loc

    def __eq__(self, other):
        return self.iNode.loc == other.iNode.loc

    def __hash__(self):
        return hash(self.iNode.loc)

    def __str__(self):
        string = f"span from {self.iNode.loc} to {self.jNode.loc}"
        return sf.fixUnits(string)
