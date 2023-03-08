# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
import collections
import operator
import sys
import numpy as np

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis.load import PointLoad
from steel_beam_analysis import units


class Node:
    """Model a node at a discrete location with associated boundary conditions,
    applied loads, demands, section properties and bracing conditions."""

    def __init__(self, loc, condition="free", **kwargs):
        # b
        self._bendingDCR = None
        self.braceBottom = kwargs.get("braceBottom", False)
        self.braceRotation = kwargs.get("braceRotation", False)
        self.braceTop = kwargs.get("braceTop", False)
        # c
        self.Cb = None
        self.condition = None
        # d
        self.__DOF = {
            "free": self.__free,
            "glass": self.__glass,
            "pin": self.__pin,
            "fix": self.__fix,
        }
        self._deflMaxAbs = {}
        self._deflMaxAbsL = {}
        # g
        self.glassDeflCheck = None
        # l
        self.leftBendingStress = 0 * units.psi
        self.leftS = None
        self.leftShearStress = 0 * units.psi
        self.loc = loc.to(units.ft)
        # m
        self.Mn = None
        self._MuMaxAbs = {}
        # p
        self.phiMn = None
        self._pointLoads = kwargs.get("pointLoads", [])
        # r
        self._rawDefls = {}
        self._rawM = {}
        self._rawMapply = {}
        self.rawMrxns = {}
        self._rawRotations = {}
        self._rawV = {}
        self._rawVapply = {}
        self.rawVrxns = {}
        self.ultVrxns = {}
        self.rotate = 0
        # s
        self._shearDCR = None
        self._span = None
        # t
        self.trans = 0
        # u
        self._ultDefls = {}
        self._ultDeflsL = {}
        self._ultM = {}
        self._ultMapply = {}
        self.ultMrxns = {}
        self._ultV = {}
        self._ultVapply = {}
        self.unbracedSpan = []
        # v
        self._VuMaxAbs = {}
        self.setCondition(condition)

    @property
    def pointLoads(self):
        return self._pointLoads

    @pointLoads.setter
    def pointLoads(self, vals):
        for pointLoad in vals:
            if not isinstance(pointLoad, PointLoad):
                raise TypeError
        self._pointLoads = vals

    @property
    def shearDCR(self):
        return self._shearDCR

    @shearDCR.setter
    def shearDCR(self, val):
        if not isinstance(val, np.float64):
            raise TypeError
        self._shearDCR = val

    @property
    def bendingDCR(self):
        return self._bendingDCR

    @bendingDCR.setter
    def bendingDCR(self, val):
        if not isinstance(val, np.float64):
            raise TypeError
        self._bendingDCR = val

    @property
    def span(self):
        return self._span

    @span.setter
    def span(self, val):
        if type(val).__name__ != "Span":
            raise TypeError
        self._span = val

    @property
    def rawMapply(self):
        return self._rawMapply

    @rawMapply.setter
    def rawMapply(self, val):
        if val is not dict:
            raise TypeError
        self._rawMapply = val

    @property
    def rawVapply(self):
        return self._rawVapply

    @rawVapply.setter
    def rawVapply(self, val):
        if val is not dict:
            raise TypeError
        self._rawVapply = val

    @property
    def rawM(self):
        return self._rawM

    @rawM.setter
    def rawM(self, val):
        if val is not dict:
            raise TypeError
        self._rawM = val

    @property
    def rawV(self):
        return self._rawV

    @rawV.setter
    def rawV(self, val):
        if val is not dict:
            raise TypeError
        self._rawV = val

    @property
    def rawDefls(self):
        return self._rawDefls

    @rawDefls.setter
    def rawDefls(self, val):
        if val and not isinstance(val, dict):
            raise TypeError
        self._rawDefls = val

    @property
    def rawRotations(self):
        return self._rawRotations

    @rawRotations.setter
    def rawRotations(self, val):
        if val is not dict:
            raise TypeError
        self._rawRotations = val

    @property
    def ultDefls(self):
        return self._ultDefls

    @ultDefls.setter
    def ultDefls(self, val):
        if val is not dict:
            raise TypeError
        self._ultDefls = val

    @property
    def deflMaxAbs(self):
        return self._deflMaxAbs

    @deflMaxAbs.setter
    def deflMaxAbs(self, val):
        if val is not dict:
            raise TypeError
        self._deflMaxAbs = val

    @property
    def ultDeflsL(self):
        return self._ultDeflsL

    @ultDeflsL.setter
    def ultDeflsL(self, val):
        if val is not dict:
            raise TypeError
        self._ultDeflsL = val

    @property
    def deflMaxAbsL(self):
        return self._deflMaxAbsL

    @deflMaxAbsL.setter
    def deflMaxAbsL(self, val):
        if val is not dict:
            raise TypeError
        self._deflMaxAbsL = val

    @property
    def ultM(self):
        return self._ultM

    @ultM.setter
    def ultM(self, val):
        if val is not dict:
            raise TypeError
        self._ultM = val

    @property
    def ultMapply(self):
        return self._ultMapply

    @ultMapply.setter
    def ultMapply(self, val):
        if val is not dict:
            raise TypeError
        self._ultMapply = val

    @property
    def MuMaxAbs(self):
        return self._MuMaxAbs

    @MuMaxAbs.setter
    def MuMaxAbs(self, val):
        if val is not dict:
            raise TypeError
        self._MuMaxAbs = val

    @property
    def ultV(self):
        return self._ultV

    @ultV.setter
    def ultV(self, val):
        if val is not dict:
            raise TypeError
        self._ultV = val

    @property
    def ultVapply(self):
        return self._ultVapply

    @ultVapply.setter
    def ultVapply(self, val):
        if val is not dict:
            raise TypeError
        self._ultVapply = val

    @property
    def VuMaxAbs(self):
        return self._VuMaxAbs

    @VuMaxAbs.setter
    def VuMaxAbs(self, val):
        if val is not dict:
            raise TypeError
        self._VuMaxAbs = val

    def addUnbracedBoundaryPt(self):
        """Take the node as an unbraced span boundary point if the compression
        flange is braced for all load combinations."""
        if self.braceTop and not self.braceBottom:
            return (
                True
                if all(val >= 0 * units.kft for val in self.ultM.values())
                else False
            )
        elif self.braceBottom and not self.braceTop:
            return (
                True
                if all(val <= 0 * units.kft for val in self.ultM.values())
                else False
            )
        if self.braceTop and self.braceBottom:
            return True
        else:
            return False

    def checkGlassDeflection(self, limit):
        """Check that the max nodal deflection doesn't exceed the limit."""
        self.glassDeflCheck = "NG" if self.deflMaxAbs["val"] >= limit else "OK"

    def assignUnbracedSpan(self, unbracedSpan):
        """Asign an unbraced span to the node."""
        self.unbracedSpan.append(unbracedSpan)

    def getMinUnbracedSpan(self):
        """If the node belongs to 2 unbraced spans, get the span with the
        minimum moment capacity. Change the data type of unbracedSpan."""

        self.unbracedSpan = min(self.unbracedSpan, key=operator.attrgetter("phiMn.val"))

    def calcReaction(self, **kwargs):
        """Calc the reaction at a node if the node is a boundary."""

        # IDEA: can what we're trying to do be done better with numpy?
        # see my stack overflow question and see if you think this is a good idea
        # https://stackoverflow.com/questions/61181688/add-items-in-numpy-arrays-where-each-item-has-an-associated-type

        leftNode = kwargs.get("leftNode", None)
        rightNode = kwargs.get("rightNode", None)
        type = kwargs.get("type", None)

        def add(*args):
            counter = collections.Counter()
            for d in args:
                counter.update(d)
            return dict(counter)

        def inv(dict):
            return {key: -val for key, val in dict.items()}

        if type == "leftEnd":
            self.rawVrxns = add(self.rawV, inv(self.rawVapply))
            self.ultVrxns = add(self.ultV, inv(self.ultVapply))
        elif type == "rightEnd":
            self.rawVrxns = add(inv(self.rawVapply), inv(self.rawV))
            self.ultVrxns = add(inv(self.ultVapply), inv(self.ultV))
        else:
            self.rawVrxns = add(rightNode.rawV, inv(leftNode.rawV), inv(self.rawVapply))
            self.ultVrxns = add(rightNode.ultV, inv(leftNode.ultV), inv(self.ultVapply))
        if self.condition == "fix":
            self.rawMrxns = add(inv(self.rawMapply), self.rawM)
            self.ultMrxns = add(inv(self.ultMapply), self.ultM)

    def setRawDispUnits(self):
        """Set units for the raw (unfactored) displacements at the node after matrix operations."""
        self.rawDefls = {k: v * units.inch for (k, v) in self.rawDefls.items()}

    def setCondition(self, condition="free"):
        """Set translational & rotational fixity for boundary conditions."""
        self.condition = self.__DOF.get(condition, self.__free)()

    def __free(self):
        self.trans = self.rotate = 0
        return "free"

    def __glass(self):
        self.trans = self.rotate = 0
        return "glass"

    def __pin(self):
        # 1 means fixed (displacement = 0), 0 means free (displacement != 0)
        self.trans = 1
        self.rotate = 0
        return "pin"

    def __fix(self):
        self.trans = self.rotate = 1
        return "fix"

    def __cmp__(self, other):
        return cmp(self.loc, other.loc)

    def __lt__(self, other):
        return self.loc < other.loc

    def __gt__(self, other):
        return self.loc > other.loc

    def __eq__(self, other):
        return round(self.loc, 5) == round(other.loc, 5)

    def __hash__(self):
        return hash(round(self.loc, 5))

    def __int__(self):
        return self.loc
