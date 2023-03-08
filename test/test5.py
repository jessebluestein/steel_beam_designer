from steel_beam_analysis.beam import *
from steel_beam_analysis import units
from steel_beam_analysis.node import *
from steel_beam_analysis.load import *
from steel_beam_analysis.steelBeam import *

node0 = Node(0 * units.ft)
node1 = Node(4.5 * units.ft, condition="pin")
node2 = Node((4.5 + 26) * units.ft, condition="pin")
node3 = Node((4.5 + 26 + 4) * units.ft)
nodes = [node0, node1, node2, node3]

iLoc = 0 * units.ft
jLoc = (4.5 + 26 + 4) * units.ft
w_D = 20 * units.psf
w_Lr = 20 * units.psf
w_S = 325 * units.psf
trib = 11 * units.ft
distLoad1 = AreaLoad(
    iLoc=iLoc, jLoc=jLoc, load=w_D, iTrib=trib, jTrib=trib, type="D", desc="roof dead"
)
distLoad2 = AreaLoad(
    iLoc=iLoc, jLoc=jLoc, load=w_Lr, iTrib=trib, jTrib=trib, type="Lr", desc="roof live"
)
distLoad3 = AreaLoad(
    iLoc=iLoc, jLoc=jLoc, load=w_S, iTrib=trib, jTrib=trib, type="S", desc="roof snow"
)
distLoads = [distLoad1, distLoad2, distLoad3]

beam = SteelBeam(
    nodes,
    rawDistLoads=distLoads,
    considerSelfWeight=True,
    depthClass=10,
    weightClass=112,
    eleSpacing=1 * units.inch,
    outputPDF=True,
    name="Typical Roof Beam",
    outputPath="output_reports/test5",
    project="Ski House",
    level="Roof",
    firm="xxx",
    engineer="xxx",
    checker="None",
)

print(beam)
