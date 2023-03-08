# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
import numpy as np
from operator import attrgetter
from sortedcontainers import SortedSet
import math

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis import units
import steel_beam_analysis.stringFixer as sf
from steel_beam_analysis.node import Node


class UnbracedSpan:
    """Model an unbraced span where the i node is the leftmost lateral bracing
    point and the j node is the rightmost lateral bracing point. The i node and
    j nodes may or may not also be vertical support points as well."""

    def __init__(self, beam, iNode, jNode):
        if (type(iNode) != Node) or (type(jNode) != Node):
            raise TypeError

        self.beam = beam
        self.cantilever = None
        self.Cb = self.Cb(self)
        self.Fcr = self.Fcr(self)
        self.iNode = iNode
        self.jNode = jNode
        self.node_qtrPt = None
        self.node_midPt = None
        self.node_threeQtrPt = None
        self._nodes = SortedSet([])
        self.len = self.jNode.loc - self.iNode.loc
        self.maxMomentNode = None
        self.phiMn = self.phiMn(self, beam)
        self.checkCantilever()

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, val):
        if not isinstance(val, Node):
            raise TypeError
        self._nodes = val

    class Cb:
        """Model lateral-torsional buckling modification factor."""

        def __init__(self, unbracedSpan):
            self.unbracedSpan = unbracedSpan
            self.ref = None
            self.val = {}

        def updateFactor(self):
            """Get/update factor and reference."""
            if self.unbracedSpan.cantilever:
                for lc in self.unbracedSpan.beam.strengthCombos:
                    self.val[str(lc)] = 1.0 * units("dimensionless")
            elif self.unbracedSpan.beam.considerCb:
                for lc in self.unbracedSpan.beam.strengthCombos:
                    MA = abs(self.unbracedSpan.node_qtrPt.ultM[str(lc)])
                    MB = abs(self.unbracedSpan.node_midPt.ultM[str(lc)])
                    MC = abs(self.unbracedSpan.node_threeQtrPt.ultM[str(lc)])
                    Mmax = abs(
                        max(
                            [node.ultM[str(lc)] for node in self.unbracedSpan.nodes],
                            key=abs,
                        )
                    )
                    if round(Mmax, 2) == 0 * units.kft:
                        self.val[str(lc)] = 1.0 * units("dimensionless")
                    else:
                        self.val[str(lc)] = (12.5 * Mmax) / (
                            2.5 * Mmax + 3 * MA + 4 * MB + 3 * MC
                        )
            else:
                for lc in self.unbracedSpan.beam.strengthCombos:
                    self.val[str(lc)] = 1.0 * units("dimensionless")
            self.ref = f"AISC/ANSI 360- {self.unbracedSpan.beam.specVersion} Sec. F1"

        def __str__(self):
            span = self.unbracedSpan
            maxCombo = self.unbracedSpan.maxMomentNode.MuMaxAbs["combo"]
            units = self.unbracedSpan.beam.outUnit["M"]
            Cb = self.val[maxCombo]
            if Cb == 1:
                return f"C_b = _bf_{{{Cb}}}"
            else:
                return sf.makeEquation(
                    "C_b",
                    "(12.5 *mul* M_max) / (2.5 *mul* M_max + 3 *mul* M_A + 4 *mul* M_B + 3 *mul* M_C)",
                    {
                        "M_{max}": round(
                            span.maxMomentNode.MuMaxAbs["val"].to(units), 1
                        ),
                        "M_A": round(span.node_qtrPt.ultM[maxCombo].to(units), 1),
                        "M_B": round(span.node_midPt.ultM[maxCombo].to(units), 1),
                        "M_C": round(span.node_threeQtrPt.ultM[maxCombo].to(units), 1),
                    },
                    round(Cb, 2),
                    lineBreak=True,
                )
                return completeEqn

    class Fcr:
        """Model the critical stress of a steel beam."""

        def __init__(self, unbracedSpan):
            self.c = 1  # for doubly-symetric I shapes
            self.ref = None
            self.unbracedSpan = unbracedSpan
            self.val = {}

        def updateFactor(self):
            """Get/update value and reference."""
            E = self.unbracedSpan.beam.E
            Lb = self.unbracedSpan.len
            self.unbracedSpan.beam.rts.updateFactor()
            rts = self.unbracedSpan.beam.rts.val
            J = self.unbracedSpan.beam.J
            c = self.c
            Sx = self.unbracedSpan.beam.Sx
            h0 = self.unbracedSpan.beam.h0
            for lc in self.unbracedSpan.beam.strengthCombos:
                Cb = self.unbracedSpan.Cb.val[str(lc)]
                self.val[str(lc)] = ((Cb * math.pi**2 * E) / (Lb / rts) ** 2) * (
                    1 + 0.078 * ((J * c) / (Sx * h0)) * (Lb / rts) ** 2
                ) ** (1 / 2)
            self.ref = f"AISC/ANSI 360- {self.unbracedSpan.beam.specVersion} Eq. F2-4"

        def __str__(self):
            maxCombo = self.unbracedSpan.maxMomentNode.MuMaxAbs["combo"]
            Cb = round(self.unbracedSpan.Cb.val[maxCombo], 2)
            variables = f"\\cfrac{{C_b \\cdot \\pi^2 \\cdot {{{{E}}}}}} {{\\left(\\cfrac{{L_b}}{{r_{{ts}}}}\\right)^2}} \\cdot \\sqrt{{1 + 0.078 \\cdot \\cfrac{{{{J}} \\cdot {{c}}}}{{{{S_x}} \\cdot h_0}} \\cdot \\left(\\cfrac{{L_b}}{{r_{{ts}}}}\\right)^2}}"
            results = variables.replace(f"C_b", str(Cb))
            results = results.replace(f"{{E}}", f"{{{self.unbracedSpan.beam.E}}}")
            results = results.replace(f"L_b", str(round(self.unbracedSpan.len, 1)))
            results = results.replace(
                f"r_{{ts}}", str(round(self.unbracedSpan.beam.rts.val, 1))
            )
            results = results.replace(f"{{J", f"{{{self.unbracedSpan.beam.J}")
            results = results.replace(f"S_x", str(self.unbracedSpan.beam.Sx))
            results = results.replace(f"h_0", str(round(self.unbracedSpan.beam.h0, 1)))
            results = results.replace(f"c}}", f"{self.c}}}")
            val = round(self.unbracedSpan.Fcr.val[maxCombo].to("ksi"), 1)
            return sf.fixUnits(
                f"F_{{cr}} & = {variables} \\\\ & = {results} = _bf_{{{val}}}"
            )

    class phiMn:
        """Model the nominal flexural strength of a steel beam."""

        def __init__(self, unbracedSpan, beam):
            self.beam = beam
            self.condition = None
            self.consideredLimitStates = None
            self.equation = None
            self.limitState = None
            self.phiB = self.phiB(self)
            self.ref = None
            self.unbracedSpan = unbracedSpan
            self.val = None
            self.variables = None

        class phiB:
            """Model the resistance factor for bending."""

            def __init__(self, phiMn):
                self.phiMn = phiMn
                self.val = None
                self.ref = None

            def updateFactor(self):
                """Get/update value and reference."""
                self.val = 0.90
                self.ref = f"AISC/ANSI 360-{self.phiMn.beam.specVersion} {{\\S}}F1a"

        def updateFactor(self):
            """Get/update value and reference."""
            self.phiB.updateFactor()
            phi = self.phiB.val
            Mp = self.beam.Mp.val
            Fy = self.beam.Fy
            Sx = self.beam.Sx
            Lb = self.unbracedSpan.len
            Lp = self.beam.Lp.val
            Lr = self.beam.Lr.val
            lam = self.beam.flange.lam.val
            lamP = self.beam.flange.lamP.val
            lamR = self.beam.flange.lamR.val
            governingLC = self.beam.maxMomentNode.MuMaxAbs["combo"]
            Fcr = self.unbracedSpan.Fcr.val[governingLC]
            Cb = self.unbracedSpan.Cb.val[governingLC]

            # limit state 0 - yielding
            ls0 = phi * Mp
            ls0eqn = f"{{\\phi_b}} \\cdot {{M_p}}"
            # limit state 1 - elastic lateral-torsional buckling (LTB)
            ls1 = phi * Cb * (Mp - (Mp - 0.7 * Fy * Sx) * ((Lb - Lp) / (Lr - Lp)))
            ls1eqn = f"{{\\phi_b}} \\cdot {{C_b}} \\cdot {{{{M_p}} - 0.7 \\cdot {{F_{{y}}}} \\cdot {{S_{{x}}}} \\cdot \\frac{{{{L_b}} - {{L_p}}}}{{{{L_r}} - {{L_p}}}}}}"
            # limit state 2 - inelastic lateral-torsional buckling (LTB)
            ls2 = phi * Fcr * Sx
            ls2eqn = f"{{\\phi_b}} \\cdot {{F_{{cr}}}} \\cdot {{S_{{x}}}}"
            # limit state 3 - flange local buckling
            ls3 = phi * (Mp - (Mp - 0.7 * Fy * Sx) * ((lam - lamP) / (lamR - lamP)))
            ls3eqn = f"{{\\phi_b}} \\cdot {{{{{{M_p}} - 0.7 \\cdot {{F_{{y}}}} \\cdot {{S_{{x}}}} \\cdot {{\\frac{{{{\\lambda}} - {{\\lambda_p}}}}{{{{\\lambda_r}} - {{\\lambda_p}}}}}}}}}}"

            if self.beam.flange.compactness == "Slender":
                sys.exit(
                    "ERROR: Input section has slender flanges and program cannot perform code checks for this condition."
                )
            if self.beam.web.compactness != "Compact":
                sys.exit(
                    "ERROR: Input section has noncompact or slender web and program cannot perform code checks for this condition."
                )

            if Lb <= Lp:
                self.condition = f"\\(_bf_{{{{L_b}} <= {{L_p}}}}\\)"
                if self.beam.flange.compactness == "Compact":
                    self.consideredLimitStates = "yielding"
                    self.equation = ls0eqn
                    self.val = ls0
                    self.limitState = "yielding"
                    self.ref = self.beam.Mp.ref
                elif self.beam.flange.compactness == "Noncompact":
                    self.consideredLimitStates = "compression flange local buckling"
                    self.equation = ls3eqn
                    self.val = ls3
                    self.limitState = "compression flange local buckling"
                    self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F3-1"
            elif Lp < Lb <= Lr:
                self.condition = f"\\(_bf_{{{{L_p}} < {{L_b}} <= {{L_r}}}}\\)"
                if self.beam.flange.compactness == "Compact":
                    self.equation = f"{ls1eqn} < {ls0eqn}"
                    self.consideredLimitStates = (
                        "LTB (not to exceed capacity based on yielding)"
                    )
                    if ls0 < ls1:
                        self.val = ls0
                        self.limitState = "yielding"
                        self.ref = self.beam.Mp.ref
                    else:
                        self.val = ls1
                        self.limitState = "LTB"
                        self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-2"
                if self.beam.flange.compactness == "Noncompact":
                    self.consideredLimitStates = "LTB (not to exceed capacity based on compression flange local buckling)"
                    self.equation = f"{ls1eqn} < {ls3eqn}"
                    if ls3 < ls1:
                        self.val = ls3
                        self.limitState = "compression flange local buckling"
                        self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F3-1"
                    else:
                        self.val = ls1
                        self.limitState = "LTB"
                        self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-2"
            elif Lb > Lr:
                self.condition = f"\\(_bf_{{{{L_b}} > {{L_r}}}}\\)"
                if self.beam.flange.compactness == "Compact":
                    self.equation = f"{ls2eqn} < {ls0eqn}"
                    self.consideredLimitStates = (
                        "inelastic LTB (not to exceed capacity based on yielding)"
                    )
                    if ls0 < ls2:
                        self.val = ls0
                        self.limitState = "yielding"
                        self.ref = self.beam.Mp.ref
                    else:
                        self.val = ls2
                        self.limitState = "inelastic LTB"
                        self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-3"
                elif self.beam.flange.compactness == "Noncompact":
                    self.consideredLimitStates = "inelsatic LTB (not to exceed capacity based on compression flange local buckling)"
                    self.equation = f"{ls2eqn} < {ls3eqn}"
                    if ls3 < ls2:
                        self.val = ls3
                        self.limitState = "compression flange local buckling"
                        self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F3-1"
                    else:
                        self.val = ls2
                        self.limitState = "inelastic LTB"
                        self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-3"

            # pass available moment to nodes in unbraced span
            for node in self.unbracedSpan.nodes:
                # for boundary nodes, phiMn is minumum considering both adjacent spans
                if node.phiMn:
                    if self.val < node.phiMn:
                        node.phiMn = self.val
                else:
                    node.phiMn = self.val

            self.condition = sf.fixUnits(self.condition)

        def __str__(self):
            span = self.unbracedSpan
            maxCombo = span.maxMomentNode.MuMaxAbs["combo"]
            results = self.equation.replace(
                f"M_p", f"{round(self.beam.Mp.val.to(units.kft), 1)}"
            )
            results = results.replace(f"F_{{y}}", str(self.beam.Fy))
            results = results.replace(f"S_{{x}}", str(round(self.beam.Sx, 1)))
            results = results.replace(f"L_p", str(round(self.beam.Lp.val, 1)))
            results = results.replace(f"L_r", str(round(self.beam.Lr.val, 1)))
            results = results.replace(
                f"F_{{cr}}", f'{round(span.Fcr.val[maxCombo].to("ksi"), 1)}'
            )
            results = results.replace(
                f"\\lambda_p", str(round(self.beam.flange.lamP.val, 1))
            )
            results = results.replace(
                f"\\lambda_r", str(round(self.beam.flange.lamR.val, 1))
            )
            results = results.replace(f"\\phi_b", str(self.phiB.val))
            results = results.replace(
                f"\\lambda", str(round(self.beam.flange.lam.val, 1))
            )
            results = results.replace(f"L_b", f"{{{round(span.len, 1)}}}")
            Cb = round(span.Cb.val[maxCombo], 2)
            results = results.replace(f"C_b", str(Cb))
            completeEqn = f"{{\\phi_b}}{{M_{{n}}}} & = {self.equation} \\\\ & = {results} = _bf_{{{round(self.beam.maxMomentNode.unbracedSpan.phiMn.val.to(units.kft), 1)}}}"
            completeEqn = sf.fixUnits(completeEqn, type="math")
            return completeEqn

    def getQuarterPointNodes(self):
        """Get the nodes closest to the quarter-point, half-point and
        three-quarter-point of the unbraced span. Nodes used for Cb calc."""

        def getActualNode(startIdx, idealPt, attrName):
            """Get the closest actual node to an ideal point on the beam."""

            diff0 = True
            diff1 = False
            while abs(diff0) > abs(diff1):
                guessLoc0 = self.beam.nodes[startIdx].loc
                diff0 = guessLoc0 - idealPt
                guessLoc1 = self.beam.nodes[int(startIdx - np.sign(diff0))].loc
                diff1 = guessLoc1 - idealPt
                startIdx -= int(np.sign(diff0))
            setattr(self, attrName, self.beam.nodes[startIdx])

        # ideal locations of three quarter points in the unbraced span
        idealQtrPt = self.iNode.loc + self.len / 4
        idealMidPt = self.iNode.loc + self.len / 2
        idealThreeQtrPt = self.iNode.loc + 3 * self.len / 4

        getActualNode(
            int(idealQtrPt / self.beam.avgNodeSpacing), idealQtrPt, "node_qtrPt"
        )
        getActualNode(
            int(idealMidPt / self.beam.avgNodeSpacing), idealMidPt, "node_midPt"
        )
        getActualNode(
            int(idealThreeQtrPt / self.beam.avgNodeSpacing),
            idealThreeQtrPt,
            "node_threeQtrPt",
        )

    def checkCantilever(self):
        """If either the i node or the j node of the span is free, report that
        the span is a cantilever."""

        supports = ["pin", "fix"]

        if (
            self.iNode.condition not in supports
            and self.iNode.loc == self.beam.nodes[0].loc
        ):
            self.cantilever = True
        elif (
            self.jNode.condition not in supports
            and self.jNode.loc == self.beam.nodes[-1].loc
        ):
            self.cantilever = True
        else:
            self.cantilever = False

    def setMaxMomentNode(self):
        """Set the node with the max absolute moment in the unbraced span."""
        self.maxMomentNode = max(self.nodes, key=lambda node: abs(node.MuMaxAbs["val"]))

    def __str__(self):
        string = f"Unbraced Span From {self.iNode.loc} to {self.jNode.loc}"
        return sf.fixUnits(string, type="text")
