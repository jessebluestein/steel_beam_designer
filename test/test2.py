from steel_beam_analysis import units
from steel_beam_analysis.load import AreaLoad, LineLoad, PointLoad
from steel_beam_analysis.node import Node
from steel_beam_analysis.steelBeam import SteelBeam
import unittest

p0 = PointLoad(shear=17.8 * units.kip, type="D", desc="2B-4 rxn")
p1 = PointLoad(shear=1.7 * units.kip, type="Lr", desc="2B-4 rxn")
p2 = PointLoad(shear=21.5 * units.kip, type="L", desc="2B-4 rxn")
p3 = PointLoad(shear=-4.6 * units.kip, type="E", desc="2B-4 rxn")
node0 = Node(9 * units.ft, pointLoads=[p0, p1, p2, p3])

p4 = PointLoad(shear=13.5 * units.kip, type="E", desc="shear wall")
node1 = Node(14.5 * units.ft, pointLoads=[p4])

node2 = Node(0 * units.ft, condition="pin")
node3 = Node(19 * units.ft, condition="pin")
node4 = Node(9.5 * units.ft)
nodes = [node0, node1, node2, node3, node4]

iLoc = 0 * units.ft
jLoc = 19 * units.ft
w0 = 23 * units.psf
w1 = 40 * units.psf
trib = 0.67 * units.ft
distLoad1 = AreaLoad(
    iLoc=iLoc, jLoc=jLoc, load=w0, iTrib=trib, jTrib=trib, type="D", desc="typ floor"
)
distLoad2 = AreaLoad(
    iLoc=iLoc, jLoc=jLoc, load=w1, iTrib=trib, jTrib=trib, type="L", desc="typ floor"
)
distLoads = [distLoad1, distLoad2]

beam = SteelBeam(
    nodes,
    rawDistLoads=distLoads,
    considerSelfWeight=True,
    depthClass=14,
    weightClass=132,
    eleSpacing=1 * units.inch,
    outputPDF=True,
    name="Test2",
    outputPath="output_reports/test2",
)

print(beam)
