from steel_beam_analysis import units
from steel_beam_analysis.load import AreaLoad, LineLoad, PointLoad
from steel_beam_analysis.node import Node
from steel_beam_analysis.steelBeam import SteelBeam

p0 = PointLoad(shear=2 * units.kip, type="D", desc="Post load")
p1 = PointLoad(shear=4 * units.kip, type="L", desc="Post load")
node0 = Node(0 * units.ft, condition="pin", pointLoads=[p0, p1])

node1 = Node(5 * units.ft, condition="glass", pointLoads=[p0, p1])
node2 = Node(15 * units.ft, condition="pin")
node3 = Node(25 * units.ft, condition="fix")
nodes = [node0, node1, node2, node3]

lineLoad0 = LineLoad(
    iLoc=0 * units.ft,
    jLoc=8 * units.yard,
    iLineLoad=2 * units.kN / units.m,
    jLineLoad=2 * units.kN / units.m,
    type="D",
    desc="Force line load (mixed units)",
)
lineLoad1 = LineLoad(
    iLoc=0 * units.ft,
    jLoc=5 * units.m,
    iLineLoad=500 * units.kg / units.m,
    jLineLoad=1000 * units.kg / units.m,
    type="L",
    desc="Mass line load (mixed units)",
)
lineLoad2 = LineLoad(
    iLoc=20 * units.ft,
    jLoc=25 * units.ft,
    iLineLoad=0.25 * units.klf,
    jLineLoad=0.25 * units.klf,
    type="L",
    desc="KLF line load",
)
lineLoads = [lineLoad0, lineLoad1, lineLoad2]

beam = SteelBeam(
    nodes,
    considerSelfWeight=False,
    rawDistLoads=lineLoads,
    depthClass=10,
    weightClass=22,
    eleSpacing=1 * units.inch,
    outputPDF=True,
    name="Test4",
    considerCb=True,
    patternLoads=["L"],
    outputPath="output_reports/test4",
)

print(beam)
