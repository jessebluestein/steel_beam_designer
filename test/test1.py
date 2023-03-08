from steel_beam_analysis.beam import *
from steel_beam_analysis import units
from steel_beam_analysis.node import *
from steel_beam_analysis.load import *
from steel_beam_analysis.steelBeam import *

p0 = PointLoad(shear=1.6 * units.kip, type="D")
p1 = PointLoad(shear=2.4 * units.kip, type="L")
node0 = Node(3 * units.ft, pointLoads=[p0, p1])

node1 = Node(0 * units.ft, condition="pin")
node2 = Node(7 * units.ft, condition="pin")
nodes = [node0, node1, node2]

beam = SteelBeam(
    nodes,
    considerSelfWeight=True,
    depthClass=12,
    weightClass=14,
    eleSpacing=0.5 * units.inch,
    outputPDF=True,
    name="Test1",
    level=2,
    outputPath="output_reports/test1",
)

print(beam)
