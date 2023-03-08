# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
from collections import defaultdict
import numpy as np

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis import units
from steel_beam_analysis.node import Node


class Element:
    """Model an element as a portion of a beam between two nodes."""

    def __init__(self, beam, iNode, jNode):
        if (type(iNode) != Node) or (type(jNode) != Node):
            raise TypeError

        self.beam = beam
        self.f0e = {}
        self._iDistLoads = {}
        self.iNode = iNode
        self._jDistLoads = {}
        self.jNode = jNode
        self.kE = None
        self.L = abs(self.jNode.loc - self.iNode.loc)
        self.formKE()

    @property
    def iDistLoads(self):
        return self._iDistLoads

    @iDistLoads.setter
    def iDistLoads(self, val):
        if val is not dict:
            sys.exit(
                f"ERROR: i-End Distributed Loads is a {type(val)} and it must be a dict."
            )
        self._iDistLoads = val

    @property
    def jDistLoads(self):
        return self._jDistLoads

    @jDistLoads.setter
    def jDistLoads(self, val):
        if val is not dict:
            sys.exit(
                f"ERROR: j-End Distributed Loads is a {type(val)} and it must be a dict."
            )
        self._jDistLoads = val

    def formf0e(self):
        """Form f0e-element body force array for trapezoidal load."""

        L = self.L.to(units.inch).magnitude

        for type in self.iDistLoads:
            w0 = self.iDistLoads[type]
            w1 = self.jDistLoads[type]

            # load array for uniform distributed load
            if w0 == w1:
                val0 = (-1 / 2) * w0 * L
                val1 = (-1 / 12) * w0 * L**2
                val2 = (-1 / 2) * w0 * L
                val3 = (-1 / 12) * w0 * L**2

            # load array for ascending trapezoidal load
            elif w0 < w1:
                val0 = (-1 / 2) * w0 * L - (3 / 20) * (w1 - w0) * L
                val1 = (-1 / 12) * w0 * L**2 - (1 / 30) * (w1 - w0) * L**2
                val2 = (-1 / 2) * w0 * L - (7 / 20) * (w1 - w0) * L
                val3 = (-1 / 12) * w0 * L**2 - (1 / 20) * (w1 - w0) * L**2

            # load array for descending trapezoidal load
            elif w0 > w1:
                val0 = (-1 / 2) * w1 * L - (7 / 20) * (w0 - w1) * L
                val1 = (-1 / 12) * w1 * L**2 - (1 / 20) * (w0 - w1) * L**2
                val2 = (-1 / 2) * w1 * L - (3 / 20) * (w0 - w1) * L
                val3 = (-1 / 12) * w1 * L**2 - (1 / 30) * (w0 - w1) * L**2

            f0e = np.array(([val0], [val1], [val2], [val3]))
            self.f0e[type] = f0e

    def formKE(self):
        """Form element stiffness array."""

        E = self.beam.E.to(units.psi).magnitude
        I = self.beam.I.to(units.inch**4).magnitude
        L = self.L.to(units.inch).magnitude

        if self.beam.considerShearDeformations:
            G = self.beam.G.to(units.psi).magnitude
            if hasattr(self.beam, "Aw"):
                # for steel beams, shear area = web area
                shearArea = self.beam.Aw.to(units.inch**2).magnitude
            else:
                # for wood beams, shear area = 5/6*cross-sectional area
                shearArea = 5 / 6 * self.beam.A.to(units.inch**2).magnitude

            b = (12 * E * I) / (G * shearArea * L**2)
            self.kE = np.matrix(
                [
                    [
                        12 * E * I / ((1 + b) * L**3),
                        6 * E * I / ((1 + b) * L**2),
                        -12 * E * I / ((1 + b) * L**3),
                        6 * E * I / ((1 + b) * L**2),
                    ],
                    [
                        6 * E * I / ((1 + b) * L**2),
                        4 * E * I * (1 + b / 4) / ((1 + b) * L),
                        -6 * E * I / ((1 + b) * L**2),
                        2 * E * I * (1 - b / 2) / ((1 + b) * L),
                    ],
                    [
                        -12 * E * I / ((1 + b) * L**3),
                        -6 * E * I / ((1 + b) * L**2),
                        12 * E * I / ((1 + b) * L**3),
                        -6 * E * I / ((1 + b) * L**2),
                    ],
                    [
                        6 * E * I / ((1 + b) * L**2),
                        2 * E * I * (1 - b / 2) / ((1 + b) * L),
                        -6 * E * I / ((1 + b) * L**2),
                        4 * E * I * (1 + b / 4) / ((1 + b) * L),
                    ],
                ]
            )

        else:
            self.kE = np.matrix(
                [
                    [
                        12 * E * I / L**3,
                        6 * E * I / L**2,
                        -12 * E * I / L**3,
                        6 * E * I / L**2,
                    ],
                    [
                        6 * E * I / L**2,
                        4 * E * I / L,
                        -6 * E * I / L**2,
                        2 * E * I / L,
                    ],
                    [
                        -12 * E * I / L**3,
                        -6 * E * I / L**2,
                        12 * E * I / L**3,
                        -6 * E * I / L**2,
                    ],
                    [
                        6 * E * I / L**2,
                        2 * E * I / L,
                        -6 * E * I / L**2,
                        4 * E * I / L,
                    ],
                ]
            )

    def __lt__(self, other):
        return self.iNode.loc < other.iNode.loc

    def __gt__(self, other):
        return self.iNode.loc > other.iNode.loc

    def __eq__(self, other):
        return self.iNode.loc == other.iNode.loc

    def __hash__(self):
        return hash(self.iNode.loc)
