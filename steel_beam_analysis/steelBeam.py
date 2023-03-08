# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
import operator, sys, os
from pytexit import py2tex
import matplotlib.pyplot as plt
import os
from datetime import date
from pdflatex import PDFLaTeX
import json

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis import units
from steel_beam_analysis.beam import *
import steel_beam_analysis.stringFixer as sf


class SteelBeam(Beam):
    """Model a steel wide flange beam."""

    def __init__(self, nodes, **kwargs):
        super().__init__(nodes, **kwargs)

        # a
        self.A = None
        self.AISC = None
        self.Aw = None
        # b
        self.bf = None
        # c
        self.compactness = None
        self.considerCb = kwargs.get("considerCb", True)
        self.considerShearDeformations = kwargs.get("considerShearDeformations", True)
        self.Cw = None
        # d
        self.d = None
        # e
        self.E = 29000 * units.ksi  # (constant for all steel)
        # f
        self.flangeCompactness = None
        self.Fu = None
        self.Fy = None
        # g
        self.G = 11200 * units.ksi  # (constant for all steel)
        self.grade = kwargs.get("grade", "A992")
        # h
        self.h = None
        self.h0 = None
        # i
        self.Ix = None
        self.Iy = None
        # j
        self.J = None
        # k
        self.kdes = None
        self.kv = self.kv(self)
        # l
        self.Lp = self.Lp(self)
        self.Lr = self.Lr(self)
        # m
        self.Mp = self.Mp(self)
        # o
        self.outputBendingCheck = None
        self.outputGeometricProps = None
        self.outputMaterialProps = None
        self.outputPath = kwargs.get("outputPath", "genericSteelBeam")
        self.outputShearCheck = None
        # p
        self.phiVn = self.phiVn(self)
        # r
        self.rts = self.rts(self)
        self.rx = None
        self.ry = None
        # s
        self.SCMversion = kwargs.get("SCMversion", 15)
        self.flange = self.Flange(self)
        self.web = self.Web(self)
        self.specVersion = None
        self.Sx = None
        self.Sy = None
        # t
        self.tf = None
        self.tw = None
        # u
        self.unitWt = 490 * units.pcf  # (constant for all steel)
        # w
        self.webCompactness = None
        self.weightClass = kwargs.get("weightClass", 100)
        # z
        self.Zx = None
        self.Zy = None

        os.makedirs(os.path.dirname(f"{self.outputPath}/test.txt"), exist_ok=True)
        self.setBeamSystem(nodes)
        self.runAnalysis()

    class Flange:
        """Model the flange of a steel beam."""

        def __init__(self, beam):
            self.beam = beam
            self.compactness = None
            self.lam = self.lamFlange(self)
            self.lamP = self.lamPflange(self)
            self.lamR = self.lamRflange(self)

        class lamFlange:
            """Model width to thickness ratio (lambda) for steel beam flange."""

            def __init__(self, flange):
                self.flange = flange
                self.ref = None
                self.val = None

            def updateFactor(self):
                b = self.flange.beam.bf / 2
                t = self.flange.beam.tf
                self.val = b / t
                self.ref = f"AISC/ANSI 360-{self.flange.beam.specVersion} Table B4.1b"
                return sf.makeEquation(
                    "\\lambda_{{flange}}",
                    "b / t",
                    {"b": round(b, 2), "t": round(t, 4)},
                    self.val,
                )

        class lamPflange:
            """Model limiting width-to-thickness ratio for flange compact/noncompact check."""

            def __init__(self, flange):
                self.flange = flange
                self.val = None
                self.ref = None

            def updateFactor(self):
                """Get/update factor and reference."""
                E = self.flange.beam.E
                Fy = self.flange.beam.Fy
                self.val = 0.38 * (E / Fy) ** (1 / 2)
                self.ref = f"AISC/ANSI 360-{self.flange.beam.specVersion} Table B4.1b"
                return (
                    f"\\lambda_{{P-flange}} = {self.flange.beam.makeE_FyString(0.38)}"
                )

        class lamRflange:
            """Model limiting width-to-thickness ratio for flange noncompact/slender check."""

            def __init__(self, flange):
                self.flange = flange
                self.val = None
                self.ref = None

            def updateFactor(self):
                """Get/update factor and reference."""
                E = self.flange.beam.E
                Fy = self.flange.beam.Fy
                self.val = 1.0 * (E / Fy) ** (1 / 2)
                self.ref = f"AISC/ANSI 360-{self.flange.beam.specVersion} Table B4.1b"
                return f"\\lambda_{{R-flange}} = {self.flange.beam.makeE_FyString(1.0)}"

        def checkCompactness(self):
            """Classify compactness of steel beam flange by comparing to limits."""
            self.lam.updateFactor()
            self.lamP.updateFactor()
            self.lamR.updateFactor()

            if self.lam.val < self.lamP.val:
                self.compactness = "Compact"
                comparison = f"\\textlambda_{{flange}} $<$ \\textlambda_{{P-flange}}"
            elif self.lamP.val < self.lam.val < self.lamR.val:
                self.compactness = "Noncompact"
                comparison = f"\\textlambda_{{P-flange}} $<$ \\textlambda_{{flange}} $<$ \\textlambda_{{R-flange}}"
            else:
                # no standard AISC shapes have slender flanges (custom sections might be allowed one day)
                self.compactness = "Slender"
                comparison = f"\\textlambda_{{R-flange}} $<$ \\textlambda_{{flange}}"

            return sf.fixUnits(
                f"{comparison} _arrow_ \; _bf_{{{self.compactness} Flange}}",
                type="text",
            )

    class kv:
        """Model web plate shear buckling coefficient."""

        def __init__(self, beam):
            self.beam = beam
            self.val = 5.34
            self.ref = None

        def updateFactor(self):
            """Get/update value and reference (value for web without transverse
            web stiffeners)."""
            self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Sec. G2.1.2"
            return sf.fixUnits(f"k_v = {self.val}")

    class Lp:
        """Model limiting laterally unbraced length for limit state of yielding."""

        def __init__(self, beam):
            self.beam = beam
            self.val = None
            self.ref = None

        def updateFactor(self):
            """Get/update value and reference."""
            self.val = (
                1.76 * self.beam.ry * (self.beam.E / self.beam.Fy) ** (1 / 2)
            ).to(units.ft)
            self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-5"
            return sf.makeEquation(
                "L_p",
                "1.76 *mul* r_y *mul* sqrt(E/F_y)",
                {"r_y": self.beam.ry, "E": self.beam.E, "F_y": self.beam.Fy},
                self.val,
            )

    class Lr:
        """Model limiting unbraced length for the limit state of inelastic LTB."""

        def __init__(self, beam):
            self.beam = beam
            self.c = 1  # for doubly-symetric I shapes
            self.val = None
            self.ref = None

        def updateFactor(self):
            """Get/update value and reference."""
            self.beam.rts.updateFactor()
            rts = self.beam.rts.val
            E = self.beam.E
            Fy = self.beam.Fy
            J = self.beam.J
            c = self.c
            Sx = self.beam.Sx
            h0 = self.beam.h0
            val0 = (((J * c) / (Sx * h0)) ** 2 + 6.76 * ((0.7 * Fy) / E) ** 2) ** (
                1 / 2
            )
            val1 = ((J * c) / (Sx * h0) + val0) ** (1 / 2)
            self.val = (1.95 * rts * (E / (0.7 * Fy)) * val1).to(units.ft)
            self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-6"
            return sf.makeEquation(
                "L_r",
                "1.95 *mul* r_ts *mul* (E/(0.7 *mul* F_y)) * sqrt((J *mul* c)/(S_x *mul* h_0) + sqrt(((J *mul* c)/(S_x *mul* h_0))**2 + 6.76 * ((0.7 *mul* F_y)/E)**2))",
                {
                    "F_y": Fy,
                    "E": E,
                    "J": round(self.beam.J, 2),
                    "c}": f"{self.c}}}",
                    "S_x": round(Sx, 1),
                    "h_0": round(h0, 1),
                    "r_{ts}": round(rts, 1),
                },
                self.val,
                lineBreak="True",
            )

    class Mp:
        """Model the plastic moment capacity of a steel beam."""

        def __init__(self, beam):
            self.beam = beam
            self.val = None
            self.ref = None

        def updateFactor(self):
            """Get/update value and reference."""
            self.val = self.beam.Fy * self.beam.Zx
            self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-1"
            return sf.makeEquation(
                "M_p",
                "F_y *mul* Z_x",
                {"F_y": self.beam.Fy, "Z_x": round(self.beam.Zx, 1)},
                self.val.to(self.beam.outUnit["M"]),
            )

    class phiVn:
        """Model the nominal shear strength of a steel beam."""

        def __init__(self, beam):
            self.beam = beam
            self.Cv1 = self.Cv1(self)
            self.val = None
            self.phiV = self.phiV(self)
            self.ref = None

        class Cv1:
            """Model the web shear strength coefficient of a steel beam."""

            def __init__(self, phiVn):
                self.phiVn = phiVn
                self.ref = None
                self.val = None

            def updateFactor(self):
                """Get/update value and reference."""
                E = self.phiVn.beam.E
                h = self.phiVn.beam.h
                tw = self.phiVn.beam.tw
                Fy = self.phiVn.beam.Fy
                self.phiVn.beam.kv.updateFactor()
                kv = self.phiVn.beam.kv.val

                if h / tw <= 2.24 * (E / Fy) ** (1 / 2):
                    comparison = "<="
                    condition = sf.makeVarString("2.24 * sqrt(E/F_y)")
                    limit = 2.24 * (E / Fy) ** (1 / 2)
                    self.val = 1.0
                    self.ref = f"AISC/ANSI 360-{self.phiVn.beam.specVersion} Eq. G2-2"
                elif h / tw <= 1.10 * ((kv * E) / Fy) ** (1 / 2):
                    comparison = "<="
                    condition = sf.makeVarString("1.10 * sqrt((k_v*E)/F_y)")
                    limit = 1.10 * ((kv * E) / Fy) ** (1 / 2)
                    self.val = 1.0
                    self.ref = f"AISC/ANSI 360-{self.phiVn.beam.specVersion} Eq. G2-3"
                elif h / tw > 1.10 * ((kv * E) / Fy) ** (1 / 2):
                    comparison = ">"
                    condition = sf.makeVarString("1.10 * sqrt((k_v*E)/F_y)")
                    limit = 1.10 * ((kv * E) / Fy) ** (1 / 2)
                    self.val = (1.10 * ((kv * E) / Fy) ** (1 / 2)) / (h / tw)
                    self.ref = f"AISC/ANSI 360-{self.phiVn.beam.specVersion} Eq. G2-4"

                vars = sf.makeVarString("h / t_w")
                ratioResults = sf.fixUnits(
                    vars, replacements={"h": round(h, 1), "t_w": round(tw, 1)}
                )
                conditionResults = sf.fixUnits(
                    condition, replacements={"E": E, "F_y": Fy, "k_v": kv}
                )
                return sf.fixUnits(
                    f"{vars} = {ratioResults} = _bf_{{{(round(h/tw, 1))}}} {comparison} {condition} = {conditionResults} = _bf_{{{round(limit, 1)}}} \\rightarrow C_{{v1}} = _bf_{{{self.val}}}"
                )

        class phiV:
            """Model the resistance factor for shear."""

            def __init__(self, phiVn):
                self.phiVn = phiVn
                self.val = None
                self.ref = None

            def updateFactor(self):
                """Get/update value and reference."""
                E = self.phiVn.beam.E
                Fy = self.phiVn.beam.Fy
                h = self.phiVn.beam.h
                tw = self.phiVn.beam.tw

                if h / tw <= 2.24 * (E / Fy) ** (1 / 2):
                    comparison = "<="
                    self.val = 1.0
                    self.ref = (
                        f"AISC/ANSI 360-{self.phiVn.beam.specVersion} {{\\S}}G2.1.a"
                    )
                else:
                    comparison = ">"
                    self.val = 0.9
                    self.ref = (
                        f"AISC/ANSI 360-{self.phiVn.beam.specVersion} {{\\S}}G1.a"
                    )

                vars = sf.makeVarString("h / t_w")
                ratioResults = sf.fixUnits(
                    vars, replacements={"h": round(h, 1), "t_w": tw}
                )
                completeEqn = f"{vars} = {ratioResults} = _bf_{{{(round(h/tw, 1))}}} {comparison} {self.phiVn.beam.makeE_FyString(2.24)} \\rightarrow \\phi_v = _bf_{{{self.val}}}"
                return sf.fixUnits(completeEqn, type="math")

        def updateFactor(self):
            """Get/update value and reference."""
            self.beam.phiVn.Cv1.updateFactor()
            self.phiV.updateFactor()
            self.val = 0.6 * self.beam.Fy * self.beam.Aw * self.beam.phiVn.Cv1.val
            self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. G2-1"
            return sf.makeEquation(
                "\\phi_v V_n",
                "0.6 *mul* F_y *mul* A_w *mul* C_v1",
                {
                    "F_y": self.beam.Fy,
                    "A_w": round(self.beam.Aw, 2),
                    "C_{v1}": round(self.beam.phiVn.Cv1.val, 2),
                    "\\phi_v": self.phiV.val,
                },
                self.val.to(self.beam.outUnit["V"]),
            )

    class rts:
        """Model rts factor used in calculation of Lr and Cb."""

        def __init__(self, beam):
            self.beam = beam
            self.val = None
            self.ref = None

        def updateFactor(self):
            """Get/update value and reference."""
            self.val = (((self.beam.Iy * self.beam.Cw) ** (1 / 2)) / self.beam.Sx) ** (
                1 / 2
            )
            self.ref = f"AISC/ANSI 360-{self.beam.specVersion} Eq. F2-7"
            return sf.makeEquation(
                "r_{{ts}}",
                "sqrt(sqrt(I_y *mul* C_w)/S_x)",
                {
                    "I_y": self.beam.Iy,
                    "C_w": self.beam.Cw,
                    "S_x": round(self.beam.Sx, 2),
                },
                self.val,
            )

    class Web:
        """Model the web of a steel beam."""

        def __init__(self, beam):
            self.beam = beam
            self.compactness = None
            self.lam = self.lamWeb(self)
            self.lamP = self.lamPweb(self)
            self.lamR = self.lamRweb(self)

        class lamWeb:
            """Model width to thickness ratio (lambda) for steel beam web."""

            def __init__(self, web):
                self.web = web
                self.ref = None
                self.val = None

            def updateFactor(self):
                self.val = self.web.beam.h / self.web.beam.tw
                self.ref = f"AISC/ANSI 360-{self.web.beam.specVersion} Table B4.1b"
                return sf.makeEquation(
                    "\\lambda_{{web}}",
                    "h/t_w",
                    {"h": round(self.web.beam.h, 2), "t_w": round(self.web.beam.tw, 3)},
                    round(self.val, 1),
                )

        class lamPweb:
            """Model limiting width-to-thickness ratio for web compact/noncompact check."""

            def __init__(self, web):
                self.web = web
                self.val = None
                self.ref = None

            def updateFactor(self):
                """Get/update factor and reference."""
                E = self.web.beam.E
                Fy = self.web.beam.Fy
                self.val = 3.76 * (E / Fy) ** (1 / 2)
                self.ref = f"AISC/ANSI 360-{self.web.beam.specVersion} Table B4.1b"
                return f"\\lambda_{{P-web}} = {self.web.beam.makeE_FyString(3.76)}"

        class lamRweb:
            """Model limiting width-to-thickness ratio for web noncompact/slender check."""

            def __init__(self, web):
                self.web = web
                self.val = None
                self.ref = None

            def updateFactor(self):
                """Get/update factor and reference."""
                E = self.web.beam.E
                Fy = self.web.beam.Fy
                self.val = 5.70 * (E / Fy) ** (1 / 2)
                self.ref = f"AISC/ANSI 360-{self.web.beam.specVersion} Table B4.1b"
                return f"\\lambda_{{R-web}} = {self.web.beam.makeE_FyString(5.70)}"

        def checkCompactness(self):
            """Classify compactness of steel beam web by comparing to limits."""
            self.lam.updateFactor()
            self.lamP.updateFactor()
            self.lamR.updateFactor()
            if self.lam.val < self.lamP.val:
                self.compactness = "Compact"
                comparison = f"\\textlambda_{{web}} $<$ \\textlambda_{{P-web}}"
            elif self.lamP.val < self.val < self.lamR.val:
                # no standard AISC shapes have noncompact web (custom sections might be allowed one day)
                self.compactness = "Noncompact"
                comparison = f"\\textlambda_{{P-web}} $<$ \\textlambda_{{web}} $<$ \\textlambda_{{R-web}}"
            else:
                # no standard AISC shapes have noncompact web (custom sections might be allowed one day)
                self.compactness = "Slender"
                comparison = f"\\textlambda_{{R-web}} $<$ \\textlambda_{{web}}"

            return sf.fixUnits(
                f"{comparison} _arrow_ \; _bf_{{{self.compactness} Web}}", type="text"
            )

    def calcCapacity(self):
        """Calculate the available flexural and shear strength (capacity)."""
        self.Mp.updateFactor()
        self.phiVn.updateFactor()
        for span in self.unbracedSpans:
            span.getQuarterPointNodes()
            span.Cb.updateFactor()
            span.Fcr.updateFactor()
            span.phiMn.updateFactor()
        for node in self.nodes:
            node.getMinUnbracedSpan()

    def calcDCRs(self):
        """Calc governing demand/capacity ratio for each node."""
        for node in self.nodes:
            node.bendingDCR = (
                (abs(node.MuMaxAbs["val"]) / node.phiMn).to_reduced_units().magnitude
            )
            node.shearDCR = (
                (abs(node.VuMaxAbs["val"]) / self.phiVn.val)
                .to_reduced_units()
                .magnitude
            )

    def checkBending(self):
        """Check that shear demand-to-capacity-ratios (DCRs) are within limits.
        Report a passing check as (OK) and a failing check as (NG)."""
        if self.maxMomentNode.bendingDCR <= 1.0:
            self.bendingCheck = "OK"
            comparison = "<"
        else:
            self.bendingCheck = "NG"
            comparison = ">"
        moment = round(abs(self.maxMomentNode.MuMaxAbs["val"].to(self.outUnit["M"])), 1)
        phiMn = round(
            self.maxMomentNode.unbracedSpan.phiMn.val.to(self.outUnit["M"]), 1
        )
        return sf.fixUnits(
            f"_bf_{{|M_u| = {moment} _sp_ {comparison} \\phi_b \\cdot M_n = {phiMn} _sp_ (DCR = {round(self.maxMomentNode.bendingDCR, 2)} - {self.bendingCheck})}}"
        )

    def checkShear(self):
        """Check that shear demand-to-capacity-ratios (DCRs) are within limits.
        Report a passing check as (OK) and a failing check as (NG)."""
        if self.maxShearNode.shearDCR <= 1.0:
            self.shearCheck = "OK"
            comparison = "<"
        else:
            self.shearCheck = "NG"
            comparison = ">"
        shear = round(abs(self.maxShearNode.VuMaxAbs["val"].to(self.outUnit["V"])), 1)
        phiVn = round(self.phiVn.val.to(self.outUnit["V"]), 1)
        return sf.fixUnits(
            f"_bf_{{|V_u| = {shear} _sp_ {comparison} \\phi_v \\cdot V_n = {phiVn} _sp_ (DCR = {round(self.maxShearNode.shearDCR, 2)} - {self.shearCheck})}}"
        )

    def checkWidthThicknessRatios(self):
        """Calculate width-to-thickness ratio and compare to limiting
        width-to-thickness ratios and determine whether web and flanges are
        compact/noncompact/slender."""
        self.flange.checkCompactness()
        self.web.checkCompactness()
        if self.flange.compactness == "Compact" and self.web.compactness == "Compact":
            self.compactness = "Compact"
        else:
            self.compactness = "Noncompact"

    def makeE_FyString(self, scalar):
        """Make a string of the square root of E/Fy multiplied by a scalar."""

        vars = sf.makeVarString(f"{scalar} *mul* sqrt(E/F_y)")
        results = sf.fixUnits(vars, replacements={"E": self.E, "F_y": self.Fy})
        val = str(
            round(
                (scalar * (self.E / self.Fy) ** (1 / 2)).to_reduced_units().magnitude, 1
            )
        )
        return f"{vars} = {results} = \\mathbf{{{val}}}"

    def outputPDFreport(self):
        """Output a PDF report to summarize analysis results."""

        # distributed loads table content
        distLoads = []
        for distLoad in self.rawDistLoads:
            startLoc = sf.fixUnits(distLoad.iLoc, type="text")
            startMag = sf.fixUnits(distLoad.iLineLoad, type="text")
            endLoc = sf.fixUnits(distLoad.jLoc, type="text")
            endMag = sf.fixUnits(distLoad.jLineLoad, type="text")
            distLoads.append(
                {
                    "startLoc": startLoc,
                    "startMag": startMag,
                    "endLoc": endLoc,
                    "endMag": endMag,
                    "type": distLoad.type,
                    "desc": distLoad.desc,
                }
            )

        # reference design vals table content
        props = [
            "Aw",
            "Cw",
            "Fu",
            "Fy",
            "Ix",
            "Iy",
            "Sx",
            "Sy",
            "Zx",
            "Zy",
            "rx",
            "ry",
            "bf",
            "tf",
            "tw",
            "h0",
            "unitWt",
        ]
        refDesignVals = {key: val for key, val in vars(self).items() if key in props}
        refDesignValsTable = []
        for prop in props:
            if len(prop) == 2:
                refDesignValsTable.append(
                    sf.fixUnits(
                        f"{prop[0]}_{{{prop[1]}}} & {round(refDesignVals[prop], 1)}",
                        type="text",
                    )
                )
        refDesignValsTable.append(
            sf.fixUnits(f"U.W. & {round(refDesignVals['unitWt'], 1)}", type="text")
        )

        # reactions table content
        reactions = ""
        for support in self.supports:
            reaction = f"{round(support.loc, 1)} & Shear & "
            for idx, type in enumerate(self.loadTypesSub):
                reaction += f'{round(support.rawVrxns[type].to(self.outUnit["V"]), 1)}'
                if idx < int(len(self.loadTypesSub) - 1):
                    reaction += " & "
                else:
                    reaction += "\\\\ \r"

            if support.condition == "fix":
                reaction += f"{round(support.loc, 1)} & Moment & "
                for idx, type in enumerate(self.loadTypesSub):
                    reaction += (
                        f'{round(support.rawMrxns[type].to(self.outUnit["M"]), 1)}'
                    )
                    if idx < int(len(self.loadTypesSub) - 1):
                        reaction += " & "
                    else:
                        reaction += "\\\\ \r"

            reaction = sf.fixUnits(reaction, type="text")
            reactions += reaction

        reactionHeaders = "Loc. & Type & "
        for idx, type in enumerate(self.loadTypesSub):
            reactionHeaders += type
            if idx < int(len(self.loadTypesSub) - 1):
                reactionHeaders += " & "
            else:
                reactionHeaders += "\\\\"

        reactionCols = ""
        for col in range(len(self.loadTypesSub) + 2):
            reactionCols += "l "

        strengthCombos = []
        for lc in self.strengthCombos:
            strengthCombos.append(
                {
                    "lc": sf.fixUnits(lc.name.replace("Lr", f"L_r"), type="text"),
                    "ref": sf.fixUnits(lc.ref, type="text"),
                }
            )

        deflCombos = []
        for lc in self.deflCombos:
            deflCombos.append(
                {
                    "lc": sf.fixUnits(lc.name.replace("Lr", f"L_r"), type="text"),
                    "ref": sf.fixUnits(lc.ref, type="text"),
                }
            )

        deflChecks = [self.maxDeflNodes["TL"].span.checkDeflections("TL")]
        if "L" in self.loadTypes:
            deflChecks.append(self.maxDeflNodes["LL"].span.checkDeflections("LL"))

        reportFile = f"{self.outputPath}/{self.projectInfo_memberName}_report.tex"
        with open(f"../steel_beam_analysis/report_template/template.tex") as f:
            template = latex_jinja_env.from_string(f.read())

        with open(f"{reportFile}", "w") as f:
            f.write(
                template.render(
                    beam=self,
                    date=date.today(),
                    distLoads=distLoads,
                    refDesignValsTable=refDesignValsTable,
                    flange=self.flange.checkCompactness(),
                    web=self.web.checkCompactness(),
                    unbracedSpans=[str(span) for span in self.unbracedSpans],
                    Lp=self.Lp.updateFactor(),
                    Lr=self.Lr.updateFactor(),
                    rts=self.rts.updateFactor(),
                    phiVn=self.phiVn.updateFactor(),
                    flangeCompactness=self.flange.compactness.lower(),
                    consideredLimitStates=self.maxMomentNode.unbracedSpan.phiMn.consideredLimitStates,
                    Mp=self.Mp.updateFactor(),
                    phiV=self.phiVn.phiV.updateFactor(),
                    Cv1=self.phiVn.Cv1.updateFactor(),
                    strengthCombos=strengthCombos,
                    deflCombos=deflCombos,
                    lamFlange=self.flange.lam.updateFactor(),
                    lamPflange=self.flange.lamP.updateFactor(),
                    lamPweb=self.web.lamP.updateFactor(),
                    lamRflange=self.flange.lamR.updateFactor(),
                    lamRweb=self.web.lamR.updateFactor(),
                    lamWeb=self.web.lam.updateFactor(),
                    spans=[str(span) for span in self.spans],
                    shearCheck=self.checkShear(),
                    bendingCheck=self.checkBending(),
                    momentsPlotFile=f"{self.outputPath}/{self.projectInfo_memberName}_moment.pgf",
                    CbPlotFile=f"{self.outputPath}/{self.projectInfo_memberName}_Cb.pgf",
                    loadDiagramPlotFile=f"{self.outputPath}/{self.projectInfo_memberName}_loadDiagram.pgf",
                    reactions=reactions,
                    reactionCols=reactionCols,
                    reactionHeaders=reactionHeaders,
                    deflChecks=deflChecks,
                    refDesignValsTableBreak1=int(len(refDesignValsTable) / 3 + 1),
                    refDesignValsTableBreak2=int(2 * len(refDesignValsTable) / 3 + 1),
                    shearsPlotFile=f"{self.outputPath}/{self.projectInfo_memberName}_shear.pgf",
                    deflsPlotFile=f"{self.outputPath}/{self.projectInfo_memberName}_defl.pgf",
                )
            )

        pdfl = PDFLaTeX.from_texfile(reportFile)
        pdfl.add_args({"-output-directory": self.outputPath})
        pdfl.create_pdf(keep_pdf_file=True, keep_log_file=False)

    def plotCb(self):
        """Plot a number line showing Cb values along the beam."""

        d = self.d.to("in").magnitude
        tf = self.tf.to("in").magnitude

        locs = [node.loc.to("ft").magnitude for node in self.unbracedSpanBoundaryPts]
        topTopFlange = [d / 2 for node in self.unbracedSpanBoundaryPts]
        botTopFlange = [d / 2 - tf for node in self.unbracedSpanBoundaryPts]
        topBotFlange = [-d / 2 for node in self.unbracedSpanBoundaryPts]
        botBotFlange = [-d / 2 + tf for node in self.unbracedSpanBoundaryPts]
        fig = plt.figure(figsize=(8, 1))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(locs, topTopFlange, color="black", linewidth=0.75)
        ax.plot(locs, botTopFlange, color="black", linewidth=0.75)
        ax.plot(locs, topBotFlange, color="black", linewidth=0.75)
        ax.plot(locs, botBotFlange, color="black", linewidth=0.75)
        for loc in locs:
            plt.axvline(x=loc, c="red", linestyle="dotted")
            # kwargs to display points: marker = 'x', markersize = 20, markeredgecolor = 'red'

        labelLocs = []
        for idx in range(len(locs) - 1):
            labelLocs.append((locs[idx + 1] + locs[idx]) / 2)

        for idx, loc in enumerate(labelLocs):
            Cb = round(
                self.unbracedSpans[idx]
                .Cb.val[self.maxMomentNode.MuMaxAbs["combo"]]
                .magnitude,
                2,
            )
            ax.text(
                loc,
                -0.025,
                f"C\\textsubscript{{b}} = {Cb}",
                bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
                ha="center",
                va="center",
            )

        plt.ylim(top=d / 2, bottom=-d / 2)

        fig.patch.set_visible(False)
        ax.axis("off")
        plt.savefig(f"{self.outputPath}/{self.projectInfo_memberName}_Cb.pgf")

    def plotMaterialSpecificFigures(self):
        """Plot figures specific to the material type."""
        self.plotCb()

    def setMaterialSpecificSystem(self):
        """Set each support as a lateral bracing point per requirements in
        AISC/ANSI 360 Sec. F1b - 'The provisions in this chapter are based on
        the assumption that points of support for beams and girders are
        restrained against rotation about their longitudinal axis'. Set
        governing unbraced length (Lp and Lr). Set resistance factor (phi)."""

        self.setLoadCombos("LRFD", self.strengthCombos)
        self.setLoadCombos("ASD", self.deflCombos)
        self.setLoadCombos("L", self.Lcombos)

        if self.bendingAxis == "strong":
            self.I = self.Ix
            self.S = self.Sx
        elif self.bendingAxis == "weak":
            self.I = self.Iy
            self.S = self.Sy
            sys.exit(
                "ERROR: Only strong axis bending currently supported for steel beams"
            )
        else:
            sys.exit(
                "ERROR: Invalid bending axis. Enter 'strong' for strong axis bending or 'weak' for weak axis bending."
            )
        self.A = self.A

        # supports restrained against rotation about beam's longitudinal axis
        for support in self.supports:
            support.braceBottom = True
            support.braceRotation = True
            support.braceTop = True
            self.unbracedSpanBoundaryPts.add(support)

        self.specVersion = 16
        self.Lp.updateFactor()
        self.Lr.updateFactor()
        self.checkWidthThicknessRatios()

    def setRefDesignValues(self):
        """Set reference design values (material and geometric properties)."""

        # get section properties from json db
        with open(f"../steel_beam_analysis/db/geometric_props.json") as f:
            geometric_props = json.load(f)
        use_section = None  # properties for the section we're using
        for section in geometric_props:
            if section["shape"] == self.shape:
                use_section = section

        # if no match, report error
        if use_section is None:
            sys.exit("ERROR: Section not found in db.")

        # get material properties from json db
        with open(f"../steel_beam_analysis/db/material_props.json") as f:
            material_props = json.load(f)
        use_material = None  # properties for the section we're using
        for material in material_props:
            if material["grade"] == self.grade:
                use_material = material

        # if no match, report error
        if use_section is None:
            sys.exit("ERROR: Material properties not found in db.")

        # define units for geometric properties
        geoPropUnits = {
            "A": units.inch**2,
            "bf": units.inch,
            "Cw": units.inch**6,
            "d": units.inch,
            "Ix": units.inch**4,
            "Iy": units.inch**4,
            "J": units.inch**4,
            "rx": units.inch,
            "ry": units.inch,
            "Sx": units.inch**3,
            "Sy": units.inch**3,
            "tf": units.inch,
            "tw": units.inch,
            "Zx": units.inch**3,
            "Zy": units.inch**3,
            "kdes": units.inch,
        }

        # define units for material properties
        matPropUnits = {"Fu": units.ksi, "Fy": units.ksi}

        # assign attributes based on queried values
        for prop, unit in geoPropUnits.items():
            self.__dict__[prop] = use_section[prop] * unit
        for prop, unit in matPropUnits.items():
            self.__dict__[prop] = use_material[prop] * unit

        self.h = self.d - 2 * self.kdes  # web height for slenderness ratio calc
        self.Aw = self.d * self.tw  # web (shear) area = OVERALL depth * web thickness
        self.h0 = self.d - self.tf  # distance between flange centroids

    def setShape(self):
        """Set shape based on depth class and weight class."""
        self.shape = f"W{self.depthClass}x{self.weightClass}"
