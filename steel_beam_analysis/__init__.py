from pint import UnitRegistry, set_application_registry
import numpy as np

units = UnitRegistry()
set_application_registry(units)
units.define("klf = kip / ft = klf")
units.define("ksf = kip / ft ** 2 = ksf")
units.define("psi = lbf / inch ** 2 = psi")
units.define("psf = lbf / ft ** 2 = psf")
units.define("pcf = lbf / ft ** 3 = pcf")
units.define("plf = lbf / ft = plf")
units.define("pli = lbf / inch = pli")
units.define("kft = kip * ft = kft")
units.define("kn = kN = kn")
units.define("lbm = lb = lbm")
units.define("lbin = lbf * inch = lbin")
units.define("jesse = 170 * lbf = jesse")


def getLoc(node, ut):
    return node.location.to(ut).magnitude


def getMag(quant):
    try:
        quant.magnitude
    except:
        return quant


def getBaseUnits(quant):
    return quant.ito_reduced_units()


vGetBaseUnits = np.vectorize(getBaseUnits)
vGetLoc = np.vectorize(getLoc)
vGetMag = np.vectorize(getMag)
