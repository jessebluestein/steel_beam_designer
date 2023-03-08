# ---------------------- OFF THE SHELF IMPORTED PACKAGES -----------------------
from pytexit import py2tex


def fixUnits(object, **kwargs):
    """Convert output pint units into a proper latex string."""

    string = str(object)
    type = kwargs.get("type", None)
    unitColor = "darkBlue"

    replacements = kwargs.get("replacements", None)
    if replacements:
        for old, new in replacements.items():
            string = string.replace(old, str(new))

    string = string.replace("kip_per_square_inch", f"{{_clr_{{_bf_{{_sp_ksi}}}}}}")
    string = string.replace("kft", f"{{_clr_{{_bf_{{_sp_kft}}}}}}")
    string = string.replace("force_pound", f"{{_clr_{{_bf_{{_sp_lb}}}}}}")
    string = string.replace("foot", f"{{_clr_{{_bf_{{_sp_ft}}}}}}")
    string = string.replace("plf", f"{{_clr_{{_bf_{{_sp_plf}}}}}}")
    string = string.replace("yard", f"{{_clr_{{_bf_{{_sp_yd}}}}}}")
    string = string.replace("kilonewton / meter", f"{{_clr_{{_bf_{{_sp_kN/m}}}}}}")
    string = string.replace(
        "kilogram * standard_gravity / meter", f"{{_clr_{{_bf_{{_sp_kg/m}}}}}}"
    )
    string = string.replace(
        "kilogram * standard_gravity", f"{{_clr_{{_bf_{{_sp_kg}}}}}}"
    )
    string = string.replace("meter", f"{{_clr_{{_bf_{{_sp_m}}}}}}")
    string = string.replace("dimensionless", "")
    string = string.replace("klf", f"{{_clr_{{_bf_{{_sp_klf}}}}}}")
    string = string.replace("kip", f"{{_clr_{{_bf_{{_sp_kip}}}}}}")
    string = string.replace("inch ** 2", f"{{_clr_{{_bf_{{_sp_inch^{{2}}}}}}}}")
    string = string.replace("inch ** 3", f"{{_clr_{{_bf_{{_sp_inch^{{3}}}}}}}}")
    string = string.replace("inch ** 4", f"{{_clr_{{_bf_{{_sp_inch^{{4}}}}}}}}")
    string = string.replace("inch ** 6", f"{{_clr_{{_bf_{{_sp_inch^{{6}}}}}}}}")
    string = string.replace("inch", f"{{_clr_{{_bf_{{_sp_in}}}}}}")
    string = string.replace("pcf", f"{{_clr_{{_bf_{{_sp_pcf}}}}}}")
    string = string.replace("_clr_", f"\\color{{{unitColor}}}")
    string = string.replace("_arrow_", "\\textrightarrow")

    if type == "text":
        string = string.replace("_bf_", "\\textbf")
        string = string.replace("_sp_", "")
        string = string.replace("^", "\\textsuperscript")
        string = string.replace("_", "\\textsubscript")
        string = string.replace("Sec. ", f"\\S")
    else:
        string = string.replace("_bf_", "\\mathbf")
        string = string.replace("_sp_", " \; ")
        string = string.replace("_clr_", f"\\color{{{unitColor}}}")
    return string


def makeVarString(string):
    """Make a latex string composed of a string of python variables."""

    string = py2tex(string, print_formula=False)
    string = string.replace("$$", "")
    string = string.replace("mul", "\cdot")
    return string


def makeEquation(name, vars, replacements, val, **kwargs):
    """Make a full equation in latex format from a python type string."""
    lineBreak = kwargs.get("lineBreak", None)
    lineBreak = "\\\\" if lineBreak else ""
    eqSign = "&=" if lineBreak else "="
    vars = makeVarString(vars)
    results = fixUnits(vars, replacements=replacements)
    val = fixUnits(str(round(val, 1)))
    return f"{name} {eqSign} {vars} {lineBreak} {eqSign} {results} {lineBreak} {eqSign} \\mathbf{{{val}}}"
