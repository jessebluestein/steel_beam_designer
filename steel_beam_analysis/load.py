# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
import sys

# -------------------------- CUSTOM IMPORTED PACKAGES --------------------------
from steel_beam_analysis import units


class LineLoad:
    def __init__(self, iLoc, jLoc, **kwargs):
        self.desc = kwargs.get("desc", "No description")
        self.iLineLoad = kwargs.get("iLineLoad", 0 * units.plf)
        self._iLoc = None
        self.iLoc = iLoc
        self.jLineLoad = kwargs.get("jLineLoad", 0 * units.plf)
        self._jLoc = None
        self.jLoc = jLoc
        self.slope = (self.jLineLoad - self.iLineLoad) / (self.jLoc - self.iLoc)
        self.type = kwargs.get("type", "D")
        self.len = self.jLoc - self.iLoc

    @property
    def iLoc(self):
        return self._iLoc

    @iLoc.setter
    def iLoc(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality != "[length]":
            sys.exit("ERROR: Line load START location must be in length units.")
        self._iLoc = val

    @property
    def jLoc(self):
        return self._jLoc

    @jLoc.setter
    def jLoc(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality != "[length]":
            sys.exit("ERROR: Line load END location must be in length units.")
        self._jLoc = val

    @property
    def iLineLoad(self):
        return self._iLineLoad

    @iLineLoad.setter
    def iLineLoad(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality == "[mass] / [time] ** 2":
            pass
        elif val.dimensionality == "[mass] / [length]":
            val *= units.gravity
        else:
            sys.exit("ERROR: Line load START loads must be in F/L or M/L units.")
        self._iLineLoad = val

    @property
    def jLineLoad(self):
        return self._jLineLoad

    @jLineLoad.setter
    def jLineLoad(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality == "[mass] / [time] ** 2":
            pass
        elif val.dimensionality == "[mass] / [length]":
            val *= units.gravity
        else:
            sys.exit("ERROR: Line load END loads must be in F/L or M/L units.")
        self._jLineLoad = val

    def __cmp__(self, other):
        return cmp(self.len, other.len)

    def __lt__(self, other):
        return self.len < other.len

    def __gt__(self, other):
        return self.len > other.len

    def __eq__(self, other):
        return self.len == other.len

    def __hash__(self):
        return hash(self.len)


class AreaLoad(LineLoad):
    """Model an area load."""

    def __init__(self, iLoc, jLoc, **kwargs):
        super().__init__(iLoc, jLoc, **kwargs)
        self._iTrib = None
        self.iTrib = kwargs.get("iTrib", 0 * units.ft)
        self._jTrib = None
        self.jTrib = kwargs.get("jTrib", 0 * units.ft)
        self._load = None
        self.load = kwargs.get("load", 0 * units.psf)
        self.iLineLoad = self.load * self.iTrib
        self.jLineLoad = self.load * self.jTrib
        self.slope = (self.jLineLoad - self.iLineLoad) / (self.jLoc - self.iLoc)

    @property
    def iTrib(self):
        return self._iTrib

    @iTrib.setter
    def iTrib(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality != "[length]":
            sys.exit("ERROR: Area load START trib must be in L units.")
        self._iTrib = val

    @property
    def jTrib(self):
        return self._jTrib

    @jTrib.setter
    def jTrib(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality != "[length]":
            sys.exit("ERROR: Area load END trib must be in L units.")
        self._jTrib = val

    @property
    def load(self):
        return self._load

    @load.setter
    def load(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality == "[mass] / [length] / [time] ** 2":
            pass
        elif valdimensionality == "[mass] / [length] ** 2":
            val *= units.gravity
        else:
            sys.exit("ERROR: Area loads must be in F/L^2 or M/L^2 units.")
        self._load = val


class PointLoad:
    """Model a point load (force or moment) applied at a node."""

    def __init__(self, type, **kwargs):
        self.desc = kwargs.get("desc", "No description")
        self._moment = None
        self.moment = kwargs.get("moment", 0 * units.kft)
        self._shear = None
        self.shear = kwargs.get("shear", 0 * units.lbf)
        self.type = type

    @property
    def shear(self):
        return self._shear

    @shear.setter
    def shear(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality == "[length] * [mass] / [time] ** 2":
            val *= -1
        elif val.dimensionality == "[mass]":
            val *= -units.gravity
        else:
            sys.exit("ERROR: All applied shears must be in force or mass units.")
        self._shear = val

    @property
    def moment(self):
        return self._moment

    @moment.setter
    def moment(self, val):
        if "Quantity" not in type(val).__name__:
            raise TypeError
        elif val.dimensionality == "[length]**2 * [mass] / [time] ** 2":
            pass
        elif val.dimensionality == "[mass] * [length]":
            val *= units.gravity
        else:
            sys.exit(
                "ERROR: All applied moments must be in force * distance or mass * distance units."
            )
        self._moment = val


class LoadCombo:
    """Model a load combination consisting of multiple unfactored loads."""

    def __init__(self, beam, loads, ref):
        self.beam = beam
        self.loads = loads
        self.ref = f"ASCE7-{self.beam.ASCE7version} {ref}"
        self.name = None
        self.makeName()
        self.setRef(ref)

    def makeName(self):
        """Make a descriptive name for the load combination."""
        name = ""
        for idx, load in enumerate(self.loads):
            name += f"{float(round(load['factor'],2))}{load['type']}"
            if idx != int(len(self.loads) - 1):
                name += " + "
        self.name = name

    def setRef(self, ref):
        """Set the reference description."""
        if ref == "L only deflection check":
            self.ref = ref
        else:
            self.ref = f"ASCE7-{self.beam.ASCE7version} {ref}"

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return f"{self.name}\n{self.ref}"
