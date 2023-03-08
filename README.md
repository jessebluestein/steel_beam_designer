# AISC Steel Beam Designer Tool

This steel beam analysis tool uses the finite element analysis (FEA) method to analyze and design wide flange beams (aka I beams or W shapes) in strong axis bending with any support conditions under any loading. The program can analyze any standard AISC shape. See [AISC Structural Steel Dimensioning Tool](https://www.aisc.org/publications/detailing-resources2/dimensioningtool/) for a list of all available W shapes. The program can handle a mix of both SI and imperial units for all input parameters. The output report is always generated in imperial units regardless of the input units.

The program generates an analysis report using [LaTeX](https://www.latex-project.org/). See the `Example Output Reports` directory for examples of the reporting format and features. The output report shows a sketch of the beam loading, along with a table of input loads. The output report also includes a bending and shear check, along with all applicable calculations in human-readable format. Finally, the program reports the deflection check and reactions at each support. Strength checks are run with ASCE 7-16 LRFD load combinations and deflection calculations are run with ASCE 7-16 ASD load combinations. For a list of load combinations and other beam properties used in design, please refer to `steel_beam_analysis/db`.

The following instructions list all of the requirements to run the program. This assumes that the user has a basic understanding of how to run a python script. [Click here](https://www.python.org/downloads/) to see the python website, where you can download the latest version of python.

# Development
From the root directory of this Github repo.

# Create virtual environment
With python installed, enter the following command into the terminal.
```
python3 -m venv venv
```

# Activate virtual environment
```
source venv/bin/activate
```

# Install requirements in virtual environment
```
pip install --upgrade pip
pip install -r requirements.txt
python setup.py install
```
Additionally, in order to render the output analysis report using LaTeX, you will need a LaTeX distribution that cannot be installed as a python package. For macOS users, you can install [MacTeX](https://www.tug.org/mactex/). Alternatively, for Windows, you should download [TeX Live Full](https://www.tug.org/texlive/acquire-netinstall.html).

# Analyze an example beam
Running the following scripts will analyze a few example beams. Also see documentation below for how to create your own beam.
```
cd test
python test1.py
python test2.py
python test3.py
python test4.py
python test5.py
```

# Creating a beam for analysis
In addition to the examples provided as test cases, please see the following instructions for how to create your own beam for analysis.

## Recognized load types
The program recognizes the following load types: `D, E, F, H, L, Lr, R, S, W`. The naming convention for load types corresponds to ASCE 7-16 variable definitions. Please refer to ASCE 7-16 for any additional descriptions of load types.

## Recognized units
The program recognizes all units in both SI and imperial systems and works with both systems interchangeably. See the following examples for how units are assigned to input parameters.

## Create Node Objects
Each beam must have nodes. Each node can be represent a pin support, fixed support, free end, or location of a load. See the next section for documentation on how to apply loads to nodes.

### Example of a free node located at 3 ft, with a couple of point loads
```
node0 = Node(3*units.ft, pointLoads = [p0, p1], desc = 'my description')
```
Note, this node is only added for the application of loads. It is not considered as a support. Loads are supplied to the node as a list. A description of the load is optional, but can be helpful for documentation in the output report.

### Example of a pinned node at 0 ft with no applied loading
```
node1 = Node (0*units.ft, condition = 'pin')
```

## Create a Point Load Object
Point load objects can be applied to beams. A couple of example point loads are as follows:

### Applied shear to beam
```
p = PointLoad(shear = 1.6 * units.kip, type = 'D')
```

### Applied moment to beam
```
p = PointLoad(moment = 3 * units.kft, type = 'D')
```

## Create an area load object
```
distLoad1 = AreaLoad(iLoc = 0 * units.ft, jLoc = 19 * units.ft, load = 23 * units.psf, iTrib = 0.67 * units.ft, jTrib = 0.67 * units.ft, type = 'D', desc = 'typ floor')
```
The area load object starts at a particular point, ends at another point, and has an associated tributary area. The tributary area can vary along the length of the load in order to create a trapezoidal loading condition.

## Create a line load object
```
lineLoad1 = LineLoad(iLoc = 0 * units.ft, jLoc = 13 * units.ft, iLineLoad = 0.026 * units.klf, jLineLoad = 0.02 * units.klf, type = 'D', desc = 'floor')
```
The line load object is similar to the area load object, except the line load's units are in force per foot, and there is no tributary area associated with the line load.

## Create a beam object

Once nodes, point loads, area loads, and/or line loads have been defined, it's time to create the beam object.

`nodes` is defined as a list of all of the nodes on the beam <br />
`rawDistLoads` is defined as a list of all line loads and area loads on the beam <br />
`considerSelfWeight` can be set to `True` or `False`. This sets whether self weight is considered as an applied load. <br />
`depthClass` is the nominal depth of the beam (i.e. for a W10x100, this would be 10) <br />
`weightClass` is the weight per foot of the beam (i.e. for a W10x100, this would be 100) <br />
`eleSpacing` is the spacing of elements in the finite element mesh. A finer mesh produces smoother moment, shear and deflection envelopes, but does not impact the overall DCRs. Using a finer mesh increases the solution time. A mesh size between 1 inch and 12 inches is appropiate for most problems. <br />
`name` is the name of the beam used in the analysis report <br />
`outputPath` is the file path where the program saves the output report <br />

An example for creating the beam is as follows:

```
beam = SteelBeam(nodes, rawDistLoads = lineLoads, considerSelfWeight = True,
depthClass = 10, weightClass = 12, eleSpacing = 2 * units.inch, outputPDF = True, name = 'Test3', outputPath = 'output_reports/test3')
```

To print out the demand capacity ratio (DCR aka stress ratio) for both bending and shear to the terminal, add the following line at the end:
```
print(beam)
```

# License

Copyright (c) 2022 Jesse Bluestein

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
