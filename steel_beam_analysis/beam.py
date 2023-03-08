# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
import numpy as np
import itertools
from sortedcontainers import SortedSet, SortedList
import warnings
import sys
import operator
import jinja2
import os
from jinja2 import Template
from pdflatex import PDFLaTeX
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import json

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis import units, vGetLoc, vGetBaseUnits, vGetMag
from steel_beam_analysis.load import AreaLoad, LineLoad, LoadCombo, PointLoad
from steel_beam_analysis.element import Element
from steel_beam_analysis.node import Node
from steel_beam_analysis.span import Span
from steel_beam_analysis.unbracedSpan import UnbracedSpan
import steel_beam_analysis.stringFixer as sf

# config jinja templating environment
latex_jinja_env = jinja2.Environment(
    block_start_string="\BLOCK{",
    block_end_string="}",
    variable_start_string="\VAR{",
    variable_end_string="}",
    comment_start_string="\#{",
    comment_end_string="}",
    line_statement_prefix="%%",
    line_comment_prefix="%#",
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader(os.path.abspath(".")),
)
latex_jinja_env.globals["len"] = len
latex_jinja_env.globals["str"] = str

# config matplotlib settings
matplotlib.use("pgf")
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
    }
)


class Beam:
    """Model a beam composed of nodes, elements, material props and loads."""

    def __init__(self, nodes, **kwargs):
        # a
        self.A = None
        self.ASCE7version = kwargs.get("ASCE7version", 16)
        self.avgNodeSpacing = None
        # b
        self.bendingAxis = kwargs.get("bendingAxis", "strong")
        self.bendingCheck = None
        # c
        self.considerSelfWeight = kwargs.get("considerSelfWeight", True)
        # d
        self.deflChecks = []
        self.deflCombos = set()
        self.deflLimGlass = kwargs.get("deflLimGlass", 0.25 * units.inch)
        self.deflRatios = kwargs.get("deflRatios", {"TL": 240, "LL": 360})
        self.depthClass = kwargs.get("depthClass", 10)
        # e
        self.elements = SortedSet([])
        self.eleSpacing = kwargs.get("eleSpacing", 1 * units.inch)
        # f
        self.F = {}
        self.F0 = {}
        self.F0Body = {}
        self.F0Node = {}
        self.FF = {}
        self.freeDOFs = []
        # g
        self.glassEverywhere = kwargs.get("glassEverywhere", False)
        # i
        self.I = None
        # k
        self.K = None
        self.KFF = None
        # l
        self.Lcombos = set()
        self.len = None
        self.loadTypes = SortedSet([])
        self.loadTypesSub = None
        # m
        self.M = {}
        self.maxMomentNode = None
        self.maxShearNode = None
        self.maxDeflNodes = {}
        # n
        self._nodes = SortedSet([])
        # o
        self.omega0 = kwargs.get("omega0", 2.5)
        self.outUnit = kwargs.get(
            "outUnit", {"M": "kft", "V": "kip", "defl": "inch", "loc": "ft"}
        )
        self.outputPDF = kwargs.get("outputPDF", False)
        self.overallCheck = None
        # p
        self.patternLoads = kwargs.get("patternLoads", ["L", "Lr"])
        self.pointLoads = []
        self.projectInfo_memberName = kwargs.get("name", "demo")
        self.projectInfo_project = kwargs.get(
            "project", "123 Maple Street, San Francisco CA"
        )
        self.projectInfo_level = kwargs.get("level", 2)
        self.projectInfo_firm = kwargs.get("firm", "ABC Company")
        self.projectInfo_engineer = kwargs.get("engineer", "Jesse")
        self.projectInfo_checker = kwargs.get("checker", "Joey")
        # r
        self._rawDistLoads = kwargs.get("rawDistLoads", [])
        self.realNodes = SortedSet([])
        self.restrainDOFs = []
        self.rho = kwargs.get("rho", 1.3)
        # s
        self.S = None
        self.SDS = kwargs.get("SDS", 1.0)
        self.seismicFactors = self.seismicFactors = {
            "omega0": self.omega0,
            "rho": self.rho,
        }
        self.seismicFactorUse = kwargs.get("seismicFactorUse", "omega0")
        self.seismicLfactor = kwargs.get("seismicLfactor", 0.5)
        self.shape = None
        self.shearCheck = None
        self.spans = SortedList()
        self.strengthCombos = set()
        self.supports = None
        # u
        self.U = {}
        self.UF = {}
        self.unbracedSpanBoundaryPts = SortedSet([])
        self.unbracedSpans = []
        # v
        self.V = {}
        # w
        self.weight = None

        # convert runtime warnings to errors to get traceback and address them
        warnings.simplefilter("error", RuntimeWarning)

    @property
    def rawDistLoads(self):
        return self._rawDistLoads

    @rawDistLoads.setter
    def rawDistLoads(self, vals):
        for distLoad in vals:
            if not isinstance(distLoad, LineLoad):
                raise TypeError
        self._rawDistLoads = vals

    def plotEnvelope(
        self,
        attr,
        units,
        combos,
        title,
        close,
        maxNodes,
        maxVals,
        labelNames,
        maxCombos,
    ):
        """Plot envelope of values at a node."""

        if close:
            locs = [0]
        else:
            locs = []

        for node in self.nodes:
            locs.append(node.loc.to("ft").magnitude)
        if close:
            locs.append(self.nodes[-1].loc.to("ft").magnitude)

        fig = plt.figure(figsize=(8, 3))
        ax = fig.add_subplot(1, 1, 1)

        for lc in combos:
            if close:
                vals = [0]
            else:
                vals = []
            for node in self.nodes:
                plotVal = getattr(node, attr)
                vals.append(plotVal[str(lc)].to(units).magnitude)
            if close:
                vals.append(0)
            ax.plot(locs, vals, linewidth=1.5)

        for idx, maxNode in enumerate(maxNodes):
            xCoord = maxNode.loc.to("ft").magnitude
            yCoord = maxVals[idx].to(units).magnitude
            maxVal = round(maxVals[idx].to(units).magnitude, 1)
            textCoord = (
                100 if xCoord <= self.len.to("ft").magnitude / 2 else -100
            )  # vary label loc
            textAlign = (
                "right" if xCoord <= self.len.to("ft").magnitude / 2 else "left"
            )  # vary label alignment
            ax.annotate(
                f"${labelNames[idx]} =$ {round(maxVals[idx].to(units), 1)}",
                xy=(xCoord, yCoord),
                xycoords="data",
                xytext=(textCoord, 0),
                textcoords="offset points",
                size=10,
                bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
                arrowprops=dict(arrowstyle="-"),
                ha=textAlign,
                va="center",
            )

            plt.text(
                0,
                0.15,
                f"Max combo: {maxCombos[idx]}",
                ha="left",
                va="top",
                transform=ax.transAxes,
                size=10,
                bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
            )

        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        plt.tick_params(
            top="off",
            bottom="off",
            left="off",
            right="off",
            labelleft="off",
            labelbottom="on",
        )
        plt.grid(b=True, which="major", color="k", linestyle="dotted")
        plt.savefig(f"{self.outputPath}/{self.projectInfo_memberName}_{title}.pgf")

    def plotLoadDiagram(self):
        """Plot a diagram of the applied loads."""

        fig, axs = plt.subplots(
            2,
            gridspec_kw={"hspace": -0.09, "height_ratios": [3, 1]},
            figsize=(8, 3),
            sharex=True,
        )
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        cycle_colors = prop_cycle.by_key()["color"]

        # create dictionaries from dist loads to store offset magnitudes
        distLoads = []
        for load in sorted(self.rawDistLoads, reverse=True):
            distLoads.append({"load": load, "offset": 0})

        # assign vertical offsets to overlapping loads
        plotLoads = []
        for load in distLoads:
            for plotLoad in plotLoads[::-1]:
                if plotLoad["load"].iLoc <= load["load"].iLoc <= plotLoad["load"].jLoc:
                    load["offset"] += (
                        max(plotLoad["load"].iLineLoad, plotLoad["load"].jLineLoad)
                        .to("plf")
                        .magnitude
                        + plotLoad["offset"]
                    )
                    break
            plotLoads.append(load)

        # plot distributed loads
        for idx, load in enumerate(plotLoads):
            iMag = load["load"].iLineLoad.to("plf").magnitude
            jMag = load["load"].jLineLoad.to("plf").magnitude
            iLoc = load["load"].iLoc.to("ft").magnitude
            jLoc = load["load"].jLoc.to("ft").magnitude
            points = [
                [iLoc, load["offset"]],
                [iLoc, iMag + load["offset"]],
                [jLoc, jMag + load["offset"]],
                [jLoc, load["offset"]],
            ]
            polygon = plt.Polygon(points, fill=True, alpha=0.4, color=cycle_colors[idx])
            axs[0].add_patch(polygon)
            axs[0].text(
                (iLoc + jLoc) / 2,
                (jMag) / 2 + load["offset"],
                f"w\\textsubscript{{{idx+1}}}",
                bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
                ha="center",
                va="center",
            )

        # plot beam flanges
        d = self.d.to("in").magnitude
        tf = d / 12  # flange width set to constant for aesthetic reasons
        locs = [0, self.len.to("ft").magnitude]
        topTopFlange = [0, 0]
        botTopFlange = [0 - tf, 0 - tf]
        topBotFlange = [-d, -d]
        botBotFlange = [-d + tf, -d + tf]
        flanges = [topTopFlange, botTopFlange, topBotFlange, botBotFlange]
        for flange in flanges:
            axs[1].plot(locs, flange, color="black", linewidth=1)
        axs[1].text(
            self.len.to("ft").magnitude / 2,
            -self.d.to("in").magnitude / 2,
            self.shape,
            bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
            ha="center",
            va="center",
        )

        # plot vertical lines at ends of beam
        leftEndX = [0, 0]
        leftEndY = [-d, 0]
        rightEndX = [self.len.to("ft").magnitude, self.len.to("ft").magnitude]
        rightEndY = leftEndY
        axs[1].plot(leftEndX, leftEndY, color="black", linewidth=1)
        axs[1].plot(rightEndX, rightEndY, color="black", linewidth=1)

        # plot gravity support locations
        pins = [support for support in self.supports if support.condition == "pin"]
        fixes = [support for support in self.supports if support.condition == "fix"]
        pinX = [pin.loc.to("ft").magnitude for pin in pins]
        pinY = [-d - 3 for pin in pins]
        fixX = [fix.loc.to("ft").magnitude for fix in fixes]
        fixY = [-d - 3 for fix in fixes]
        axs[1].scatter(pinX, pinY, marker="^", s=200, c="red")
        axs[1].scatter(fixX, fixY, marker="s", s=200, c="blue")

        # plot dimensions between supports
        spanPts = [span.iNode.loc.to("ft").magnitude for span in self.spans]
        spanPts.append(self.len.to("ft").magnitude)
        for idx, support in enumerate(spanPts):
            if idx != len(spanPts) - 1:
                dist = spanPts[idx + 1] - support
                # plot dimension line (no text)
                axs[1].annotate(
                    f"",
                    xy=(support, -d - 5),
                    xycoords="data",
                    xytext=(support + dist, -d - 5),
                    textcoords="data",
                    arrowprops=dict(arrowstyle="<->, head_width=0.5", color="#33ADFF"),
                    ha="center",
                )
                # plot text in center of dimension line
                axs[1].text(
                    support + dist / 2,
                    -d - 5,
                    f"Span {idx} = {dist} ft",
                    bbox=dict(boxstyle="round4, pad=0.5", fc="0.8"),
                    size=10,
                    ha="center",
                    va="center",
                )

        # plot applied point loads
        pointLoadLocs = []
        pointLoadNodes = [node for node in self.nodes if node.pointLoads]
        for node in pointLoadNodes:
            for pointLoad in node.pointLoads:
                pointLoadLocs.append(node.loc.to("ft").magnitude)
        pointLoadLocs = set(pointLoadLocs)
        for loc in pointLoadLocs:
            axs[0].annotate(
                f"$P_x$",
                xy=(loc, 0),
                xycoords="data",
                xytext=(0, 100),
                textcoords="offset points",
                size=12,
                bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
                arrowprops=dict(arrowstyle="->, head_width=0.5"),
                ha="center",
            )

        # plot settings and save
        fig.patch.set_visible(False)
        axs[0].axis("off")
        axs[1].axis("off")
        axs[1].autoscale()
        axs[0].autoscale()
        plt.savefig(
            f"{self.outputPath}/{self.projectInfo_memberName}_loadDiagram.pgf",
            dpi=90,
            pad_inches=0.5,
        )

    def runAnalysis(self):
        """Run analysis on the beam system."""

        # ---------------------- FORM ELEMENT LOAD ARRAYS ----------------------
        for element in self.elements:
            element.formf0e()

        # ------------------------ CALC DEMAND USING FEA -----------------------
        # demands from applied distributed loads on elements
        for type in self.loadTypesSub:
            self.F0Body[type] = np.zeros((len(self.nodes) * 2, 1))
            for idx, element in enumerate(self.elements):
                f0e = element.f0e.get(type, np.array(([0], [0], [0], [0])))
                self.F0Body[type][2 * idx : 2 * idx + 4] = np.add(
                    self.F0Body[type][2 * idx : 2 * idx + 4], f0e
                )

        # form F0Node
        for type in self.loadTypesSub:
            self.F0Node[type] = np.zeros((len(self.nodes) * 2, 1))
            for idx, node in enumerate(self.nodes):
                self.F0Node[type][2 * idx, 0] = (
                    node.rawVapply.get(type, 0 * units.lbf).to(units.lbf).magnitude
                )
                self.F0Node[type][2 * idx + 1, 0] = (
                    node.rawMapply.get(type, 0 * units.lbin).to(units.lbin).magnitude
                )

        # form point loads list used for plotting and output table
        idx = 1  # starts at 1 because used in output
        for node in self.nodes:
            for pointLoad in node.pointLoads:
                loc = sf.fixUnits(node.loc, type="text")
                shear = sf.fixUnits(-pointLoad.shear, type="text")
                self.pointLoads.append(
                    {
                        "id": f"P\\textsubscript{{{idx}}}",
                        "loc": loc,
                        "shear": shear,
                        "type": pointLoad.type,
                        "desc": pointLoad.desc,
                    }
                )
                idx += 1

        # combination of demands from elements and nodes
        for type in self.loadTypesSub:
            self.F0[type] = np.add(self.F0Body[type], self.F0Node[type])

        # global applied forces at free DOFs
        for type in self.loadTypesSub:
            self.FF[type] = self.F0[type][np.ix_(self.freeDOFs)]

        # global stiffness array
        self.K = np.zeros((len(self.nodes) * 2, len(self.nodes) * 2))
        for idx, element in enumerate(self.elements):
            self.K[2 * idx : 2 * idx + 4, 2 * idx : 2 * idx + 4] = np.add(
                self.K[2 * idx : 2 * idx + 4, 2 * idx : 2 * idx + 4], element.kE
            )

        # global stiffness at free DOFs
        self.KFF = self.K[np.ix_(self.freeDOFs, self.freeDOFs)]

        # displacements at free DOFs
        for type in self.loadTypesSub:
            self.UF[type] = np.matmul(np.linalg.inv(self.KFF), self.FF[type])

        # pass displacements at free DOFs to nodes
        for type in self.loadTypesSub:
            for idx, dof in enumerate(self.freeDOFs):
                if dof % 2 == 0:
                    self.nodes[int(dof / 2)].rawDefls[type] = self.UF[type][idx, 0]
                else:
                    self.nodes[int((dof - 1) / 2)].rawRotations[type] = self.UF[type][
                        idx, 0
                    ]

        # assemble displacement array (translational + rotational) for global system
        for type in self.loadTypesSub:
            self.U[type] = np.zeros((len(self.nodes) * 2, 1))
            for idx, node in enumerate(self.nodes):
                self.U[type][2 * idx, 0] = node.rawDefls.get(type, 0)
                self.U[type][2 * idx + 1, 0] = node.rawRotations.get(type, 0)

        # set units for raw displacements at nodes
        for node in self.nodes:
            node.setRawDispUnits()

        # forces (moments + shears) in global system
        for type in self.loadTypesSub:
            F = np.zeros((len(self.nodes) * 2, 1))
            fEvals = []
            for idx, elem in enumerate(self.elements):
                f0e = elem.f0e.get(type, np.array(([0], [0], [0], [0])))
                fEvals.append(
                    np.add(
                        np.matmul(elem.kE, self.U[type][2 * idx : 2 * idx + 4]), -f0e
                    )
                )
            for idx, elem in enumerate(self.elements):
                F[2 * idx : 2 * idx + 2] = fEvals[idx][0:2]
            F[-2:] = -fEvals[len(self.elements) - 1][2:4]
            self.F[type] = F

        # extract bending and shear demands from global array, set nodal values
        for type in self.loadTypesSub:
            self.M[type] = -self.F[type][1::2]  # sign flipped for convention
            self.V[type] = self.F[type][0::2]
            for idx, node in enumerate(self.nodes):
                node.rawM[type] = self.M[type][idx, 0] * units.lbin
                node.rawV[type] = self.V[type][idx, 0] * units.lbf

        # calc factored values and set values at each node
        # IDEA: any "shorthand" way to do all these nested for loops?
        # https://stackoverflow.com/questions/5462047/pythonic-shortcut-for-doubly-nested-for-loops
        for node in self.nodes:
            for lc in self.strengthCombos:
                node.ultV[str(lc)] = sum(
                    [
                        node.rawV.get(load["type"], 0) * load["factor"]
                        for load in lc.loads
                    ]
                )
                node.ultVapply[str(lc)] = sum(
                    [
                        node.rawVapply.get(load["type"], 0) * load["factor"]
                        for load in lc.loads
                    ]
                )
                node.ultM[str(lc)] = sum(
                    [
                        node.rawM.get(load["type"], 0) * load["factor"]
                        for load in lc.loads
                    ]
                )
                node.ultMapply[str(lc)] = sum(
                    [
                        node.rawMapply.get(load["type"], 0) * load["factor"]
                        for load in lc.loads
                    ]
                )
            for lc in self.deflCombos:
                node.ultDefls[str(lc)] = sum(
                    [
                        node.rawDefls.get(load["type"], 0 * units.inch) * load["factor"]
                        for load in lc.loads
                    ]
                )
            for lc in self.Lcombos:
                node.ultDeflsL[str(lc)] = sum(
                    [
                        node.rawDefls.get(load["type"], 0 * units.inch) * load["factor"]
                        for load in lc.loads
                    ]
                )
            node.deflMaxAbs["combo"] = max(
                node.ultDefls, key=lambda y: abs(node.ultDefls[y])
            )
            node.deflMaxAbs["val"] = node.ultDefls[node.deflMaxAbs["combo"]]
            node.MuMaxAbs["combo"] = max(node.ultM, key=lambda y: abs(node.ultM[y]))
            node.MuMaxAbs["val"] = node.ultM[node.MuMaxAbs["combo"]]
            node.VuMaxAbs["combo"] = max(node.ultV, key=lambda y: abs(node.ultV[y]))
            node.VuMaxAbs["val"] = node.ultV[node.VuMaxAbs["combo"]]

            if "L" in self.loadTypes:
                node.deflMaxAbsL["combo"] = max(
                    node.ultDeflsL, key=lambda y: abs(node.ultDeflsL[y])
                )
                node.deflMaxAbsL["val"] = node.ultDeflsL[node.deflMaxAbsL["combo"]]

        # get max demand nodes
        self.maxMomentNode = max(self.nodes, key=lambda node: abs(node.MuMaxAbs["val"]))
        self.maxShearNode = max(self.nodes, key=lambda node: abs(node.VuMaxAbs["val"]))
        self.maxDeflNodes["TL"] = max(
            self.nodes, key=lambda node: abs(node.deflMaxAbs["val"])
        )
        if "L" in self.loadTypes:
            self.maxDeflNodes["LL"] = max(
                self.nodes, key=lambda node: abs(node.deflMaxAbsL["val"])
            )

        # set max total load & live load deflection nodes in each span
        for span in self.spans:
            span.setMaxTLdeflNode()
            if "L" in self.loadTypes:
                span.setMaxLLdeflNode()

        # -------------------------- SET UNBRACED SPANS ------------------------
        # always include both ends of the beam (even if cantilevered ends)
        self.unbracedSpanBoundaryPts.add(self.nodes[0])
        self.unbracedSpanBoundaryPts.add(self.nodes[-1])
        # add in top/bottom brace points based on sign of moment demand at node
        for node in self.nodes:
            if node.addUnbracedBoundaryPt():
                self.unbracedSpanBoundaryPts.add(node)
        for idx, node in enumerate(self.unbracedSpanBoundaryPts):
            if idx != 0:
                self.unbracedSpans.append(
                    UnbracedSpan(self, self.unbracedSpanBoundaryPts[idx - 1], node)
                )
        spanIter = 0
        for node in self.nodes:
            if node.loc == self.unbracedSpans[spanIter].jNode.loc:
                self.unbracedSpans[spanIter].nodes.add(node)
                node.assignUnbracedSpan(self.unbracedSpans[spanIter])
                spanIter += 1
            if spanIter < len(self.unbracedSpans):
                self.unbracedSpans[spanIter].nodes.add(node)
                node.assignUnbracedSpan(self.unbracedSpans[spanIter])

        # set max moment in unbraced spans
        for span in self.unbracedSpans:
            span.setMaxMomentNode()

        # ------------- CALC CAPACITY, DCRs AND CHECK BENDING/SHEAR ------------
        self.calcCapacity()
        self.calcDCRs()
        self.checkBending()
        self.checkShear()

        # --------------------------- CALC REACTIONS ---------------------------
        for idx, node in enumerate(self.nodes):
            if node in self.supports:
                if node.loc == 0 * units.ft:
                    node.calcReaction(type="leftEnd")
                elif node.loc == self.len:
                    node.calcReaction(type="rightEnd")
                else:
                    node.calcReaction(
                        leftNode=self.nodes[idx - 1], rightNode=self.nodes[idx + 1]
                    )

        # ------------------------- CHECK DEFLECTIONS --------------------------
        for span in self.spans:
            if span.maxDeflNodes["TL"].deflMaxAbs["val"] != 0 * units.inch:
                span.setTLdeflRatios(self.deflRatios["TL"])
            span.checkDeflections("TL")
            if "L" in self.loadTypes:
                if span.maxDeflNodes["LL"].deflMaxAbsL["val"] != 0 * units.inch:
                    span.setLLdeflRatios(self.deflRatios["LL"])
                    span.checkDeflections("LL")

        for span in self.spans:
            self.deflChecks.append(span.isDeflectionOK())

        # deflection check w/ glass at discrete points
        for node in self.nodes:
            if node.condition == "glass":
                node.checkGlassDeflection(self.deflLimGlass)

        # set beam deflection check at glass based on nodal deflections
        self.deflChecks.append(
            "NG" if any(n.glassDeflCheck == "NG" for n in self.nodes) else "OK"
        )

        # if glass everywhere, check that glass deflection passes everywhere
        if self.glassEverywhere:
            self.deflChecks.append(
                "OK"
                if all(abs(n.deflMaxAbs["val"]) < self.deflLimGlass for n in self.nodes)
                else "NG"
            )

        # --------------------------- CHECK OVERALL ----------------------------
        checks = []
        checks.append(self.bendingCheck)
        checks.append(self.shearCheck)
        for check in self.deflChecks:
            checks.append(check)

        self.overallCheck = "OK" if all(c == "OK" for c in checks) else "NG"

        # ----------------------------- OUTPUT PDF  ----------------------------
        maxDeflNodes = [self.maxDeflNodes["TL"]]
        maxDeflVals = [self.maxDeflNodes["TL"].deflMaxAbs["val"]]
        maxDeflLabels = ["\Delta"]
        maxDeflAnnos = [self.maxDeflNodes["TL"].deflMaxAbs["combo"]]
        if "L" in self.loadTypes:
            maxDeflNodes.append(self.maxDeflNodes["LL"])
            maxDeflVals.append(self.maxDeflNodes["LL"].deflMaxAbsL["val"])
            maxDeflLabels.append("\Delta")
            maxDeflAnnos.append(self.maxDeflNodes["LL"].deflMaxAbsL["combo"])

        if self.outputPDF:
            self.plotLoadDiagram()

            self.plotEnvelope(
                "ultV",
                self.outUnit["V"],
                self.strengthCombos,
                "shear",
                True,
                [self.maxShearNode],
                [self.maxShearNode.VuMaxAbs["val"]],
                ["V_u"],
                [self.maxShearNode.VuMaxAbs["combo"]],
            )

            self.plotEnvelope(
                "ultM",
                self.outUnit["M"],
                self.strengthCombos,
                "moment",
                True,
                [self.maxMomentNode],
                [self.maxMomentNode.MuMaxAbs["val"]],
                ["M_u"],
                [self.maxMomentNode.MuMaxAbs["combo"]],
            )

            self.plotEnvelope(
                "ultDefls",
                self.outUnit["defl"],
                self.deflCombos,
                "defl",
                False,
                maxDeflNodes,
                maxDeflVals,
                maxDeflLabels,
                maxDeflAnnos,
            )

            self.plotMaterialSpecificFigures()
            self.outputPDFreport()

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, vals):
        for node in vals:
            if not isinstance(node, Node):
                raise TypeError
        self._nodes = vals

    def setBeamSystem(self, nodes):
        """Set generic beam system parameters."""

        if len(set(nodes)) < len(nodes):
            sys.exit("ERROR: Multiple nodes cannot have the same location.")

        for node in nodes:
            self.nodes.add(node)

        self.setShape()
        self.setRefDesignValues()
        self.len = self.nodes[-1].loc - self.nodes[0].loc

        # ----------------- SET SPAN GEOMETRY (i AND j NODES) ------------------
        for idx, node in enumerate(self.nodes):
            if idx == 0:
                self.realNodes.add(node)  # 'real' nodes define spans
            elif node.condition == "pin" or node.condition == "fix":
                self.realNodes.add(node)
            elif idx == (len(self.nodes) - 1):
                self.realNodes.add(node)

        for idx, iNode in enumerate(self.realNodes[:-1]):
            self.spans.add(Span(iNode, self.realNodes[idx + 1]))

        supportTypes = ["pin", "fix"]
        self.supports = [
            node for node in self.realNodes if node.condition in supportTypes
        ]

        # --------------------------- SET LOAD TYPES ---------------------------
        # NOTE: if load types are set after distributed loads, then seperate
        # lines for including D if self weight is to be considered can be removed

        # include self weight in load types if self weight is considered
        if self.considerSelfWeight == True:
            self.loadTypes.add("D")

        # include load types for point loads
        for node in self.nodes:
            for pointLoad in node.pointLoads:
                if not isinstance(pointLoad, PointLoad):
                    raise TypeError
                self.loadTypes.add(pointLoad.type)

        # include load types for distributed loads
        for distLoad in self.rawDistLoads:
            self.loadTypes.add(distLoad.type)

        # remove pattern load types that aren't actually on beam
        self.patternLoads = [
            load for load in self.patternLoads if load in self.loadTypes
        ]

        # set subdivided load types list (non-pattern load types + pattern load types w/ indicies)
        tempList0 = [load for load in self.loadTypes if load not in self.patternLoads]
        tempList1 = []
        for patLoad in self.patternLoads:
            for idx, span in enumerate(self.spans):
                tempList1.append(f"{patLoad}{idx}")
        self.loadTypesSub = tempList0 + tempList1

        self.setMaterialSpecificSystem()

        # ----------- CREATE NODES THAT ARE PROGRAMATICALLY REQUIRED -----------
        # create nodes where dist loads start and end if no node already exists
        if self.rawDistLoads:
            for distLoad in self.rawDistLoads:
                self.nodes.add(Node(distLoad.iLoc))
                self.nodes.add(Node(distLoad.jLoc))

        # add mesh nodes at element spacing interval
        meshNodes = []
        for idx in range(len(self.nodes) - 1):
            # skip the nodes it pertains to
            for location in np.arange(
                self.nodes[idx].loc + self.eleSpacing,
                self.nodes[idx + 1].loc - self.eleSpacing,
                self.eleSpacing,
                dtype="object",
            ):
                meshNodes.append(Node(location))
        for meshNode in meshNodes:
            self.nodes.add(meshNode)

        # ------------------------ ASSIGN NODES TO SPANS -----------------------
        spanIter = 0
        for node in self.nodes:
            if node.loc == self.spans[spanIter].jNode.loc:
                self.spans[spanIter].nodes.add(node)
                spanIter += 1
            if spanIter < len(self.spans):
                self.spans[spanIter].nodes.add(node)
                node.span = self.spans[spanIter]

        # ----------------------- CREATE ELEMENT GEOMETRY ----------------------
        for idx, iNode in enumerate(self.nodes[:-1]):
            self.elements.add(Element(self, iNode, self.nodes[idx + 1]))

        # ----------------------- ASSIGN ELEMENTS TO SPANS ---------------------
        spanIter = 0
        for element in self.elements:
            if element.jNode.loc == self.spans[spanIter].jNode.loc:
                self.spans[spanIter].elements.add(element)
                spanIter += 1
            if spanIter < len(self.spans):
                self.spans[spanIter].elements.add(element)

        # ----------- CALC SELF WEIGHT & INCLUDE AS DISTRIBUTED LOAD -----------
        self.weight = round((self.unitWt * self.A).to(units.plf), 1)
        if self.considerSelfWeight:
            self.rawDistLoads.append(
                LineLoad(
                    iLoc=0 * units.ft,
                    jLoc=self.len,
                    iLineLoad=self.weight,
                    jLineLoad=self.weight,
                    desc="Self weight",
                )
            )

        # ----------------- ASSIGN SUPERIMPOSED LOADS TO ELEMENTS --------------
        if self.rawDistLoads:
            for dl in self.rawDistLoads:
                rangeElems = [
                    elem
                    for elem in self.elements
                    if elem.iNode.loc >= dl.iLoc and elem.jNode.loc <= dl.jLoc
                ]
                for elem in rangeElems:
                    iDist = elem.iNode.loc - dl.iLoc
                    jDist = elem.jNode.loc - dl.iLoc
                    iMag = dl.iLineLoad + dl.slope * (iDist)
                    jMag = dl.iLineLoad + dl.slope * (jDist)
                    elem.iDistLoads[dl.type] = (
                        elem.iDistLoads.get(dl.type, 0) + iMag.to(units.pli).magnitude
                    )
                    elem.jDistLoads[dl.type] = (
                        elem.jDistLoads.get(dl.type, 0) + jMag.to(units.pli).magnitude
                    )

        # ------------------- ASSIGN SUPERIMPOSED LOADS TO NODES ---------------
        for node in self.nodes:
            for type in [pointLoad.type for pointLoad in node.pointLoads]:
                node.rawVapply[type] = sum(
                    [ptLd.shear for ptLd in node.pointLoads if ptLd.type == type]
                )
                node.rawMapply[type] = sum(
                    [ptLd.moment for ptLd in node.pointLoads if ptLd.type == type]
                )

        # ------------------------ BREAK OUT PATTERN LOADS ---------------------
        # for applied distributed loads
        for idx, span in enumerate(self.spans):
            for elem in span.elements:
                for patLoad in self.patternLoads:
                    elem.iDistLoads[f"{patLoad}{idx}"] = elem.iDistLoads.pop(patLoad, 0)
                    elem.jDistLoads[f"{patLoad}{idx}"] = elem.jDistLoads.pop(patLoad, 0)

        # for applied point loads
        for idx, span in enumerate(self.spans):
            for node in span.nodes:
                for patLoad in self.patternLoads:
                    node.rawVapply[f"{patLoad}{idx}"] = node.rawVapply.pop(
                        patLoad, 0 * units.kip
                    )
                    node.rawMapply[f"{patLoad}{idx}"] = node.rawMapply.pop(
                        patLoad, 0 * units.kft
                    )

        # ------- SET SYSTEM PARAMETERS AND CHECK THAT ANALYSIS CAN RUN --------
        for idx, node in enumerate(self.nodes):
            if node.trans:
                self.restrainDOFs.append(2 * idx)
            else:
                self.freeDOFs.append(2 * idx)
            if node.rotate:
                self.restrainDOFs.append(2 * idx + 1)
            else:
                self.freeDOFs.append(2 * idx + 1)

        self.avgNodeSpacing = self.len / len(self.nodes)

        if len(self.restrainDOFs) <= 1:
            sys.exit(
                "ERROR: Insufficient supports provided. Beam needs (2) pinned nodes or (1) fixed node to be stable."
            )

        for idx, node in enumerate(self.nodes):
            if node.condition == "fix":
                if idx == 0 or idx == len(self.nodes) - 1:
                    continue
                else:
                    sys.exit("ERROR: Fixed nodes can only be at beginning or end.")

    def setLoadCombos(self, collection, targetAttr):
        """Set load combinations for strength and deflection checks given a
        collection of load combinations to pull from and a target attribute to
        set with the list of load combinations."""

        # determine which collection to look for load combos
        if collection == "LRFD":
            db_path = f"../steel_beam_analysis/db/lrfd_combos.json"
        elif collection == "ASD":
            db_path = f"../steel_beam_analysis/db/asd_combos.json"
        elif collection == "L":
            db_path = f"../steel_beam_analysis/db/L_combos.json"
        else:
            sys.exit("bad collection option for load combos!")

        # read load combo data from json db
        with open(db_path) as f:
            raw_combos = json.load(f)

        # filter raw load combo data
        filtered_combos = []
        for combo in raw_combos:
            filtered_combo = {}
            for k in combo.keys():
                if combo[k] != None and k in self.loadTypes:
                    filtered_combo[k] = combo[k]
                if k == "ref":
                    filtered_combo[k] = combo[k]
            if len(filtered_combo) > 1:
                filtered_combos.append(filtered_combo)  # don't append 0 length combos

        # build load combo objects
        for combo in filtered_combos:
            comboNoRef = {k: v for k, v in combo.items() if k != "ref"}
            patLoads = list(set(list(comboNoRef.keys())) & set(self.patternLoads))
            nonPatLoads = [load for load in comboNoRef if load not in patLoads]
            if patLoads:
                tfPerms = list(itertools.product([True, False], repeat=len(self.spans)))
                spanIdxsCombos = []
                for perm in tfPerms:
                    spanIdxsCombos.append([i for i, v in enumerate(perm) if v])
                spanIdxsCombos.remove([])
                for perm in spanIdxsCombos:
                    lcLoads = []
                    for load in nonPatLoads:
                        lcLoads.append({"type": load, "factor": eval(str(combo[load]))})
                    for spanIdx in perm:
                        for load in patLoads:
                            lcLoads.append(
                                {
                                    "type": f"{load}{spanIdx}",
                                    "factor": eval(str(combo[load])),
                                }
                            )
                    if lcLoads:
                        targetAttr.add(LoadCombo(self, lcLoads, combo["ref"]))
            else:
                lcLoads = []
                for load in nonPatLoads:
                    lcLoads.append({"type": load, "factor": eval(str(combo[load]))})
                if lcLoads:
                    targetAttr.add(LoadCombo(self, lcLoads, combo["ref"]))

    def __str__(self):
        bendingDCR = round(self.maxMomentNode.bendingDCR, 3)
        shearDCR = round(self.maxShearNode.shearDCR, 3)
        string = f"Bending DCR... \t{bendingDCR}\n"
        string += f"Shear DCR... \t{shearDCR}\n"
        string += f"Analysis ran successfully!"
        return string
