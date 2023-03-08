from steel_beam_analysis import units
from steel_beam_analysis.load import AreaLoad, LineLoad, PointLoad
from steel_beam_analysis.node import Node
from steel_beam_analysis.steelBeam import SteelBeam

p0 = PointLoad(shear=0.39 * units.kip, type="D")
p1 = PointLoad(shear=0.34 * units.kip, type="Lr")
p2 = PointLoad(shear=0.01 * units.kip, type="L")
node0 = Node(0 * units.ft, condition="free", pointLoads=[p0, p1, p2])

p2 = PointLoad(shear=1.3 * units.kip, type="D")
p3 = PointLoad(shear=0.7 * units.kip, type="Lr")
p4 = PointLoad(shear=0.3 * units.kip, type="L")
node1 = Node(6 * units.ft, condition="pin", pointLoads=[p2, p3, p4])

p5 = PointLoad(shear=3.57 * units.kip, type="D")
p6 = PointLoad(shear=0.64 * units.kip, type="Lr")
p7 = PointLoad(shear=2.67 * units.kip, type="L")
p8 = PointLoad(shear=0.11 * units.kip, type="E")
node2 = Node(13 * units.ft, condition="pin", pointLoads=[p5, p6, p7, p8])
nodes = [node0, node1, node2]

w0mag = 0.026 * units.klf
w1mag = 0.02 * units.klf
iLoc = 0 * units.ft
jLoc = 13 * units.ft
lineLoad1 = LineLoad(
    iLoc=iLoc, jLoc=jLoc, iLineLoad=w0mag, jLineLoad=w0mag, type="D", desc="floor"
)
lineLoad2 = LineLoad(
    iLoc=iLoc, jLoc=jLoc, iLineLoad=w1mag, jLineLoad=w1mag, type="Lr", desc="floor"
)
lineLoads = [lineLoad1, lineLoad2]

beam = SteelBeam(
    nodes,
    rawDistLoads=lineLoads,
    considerSelfWeight=True,
    depthClass=10,
    weightClass=12,
    eleSpacing=2 * units.inch,
    bendingAxis="strong",
    outputPDF=True,
    name="Test3",
    outputPath="output_reports/test3",
)

print(beam)
