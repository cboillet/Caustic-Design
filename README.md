# Caustic-Design

### 1. Image->Singularities 
// Assigned to: Cam<br>
<b>Input</b>
 *  Image (PNG)

<b>Output</b>
 *  Singularities (SVG-File) [1]

### 2. Irradiance-> Vornoi diagramm Python or CGAL 
// Assigned to: Patrick<br>
<b>Input</b>
 *  Singularities [1]
 *  Area Distribution

<b>Output</b>
 *  Voronoi Diagram [2]<br>

<b>Building Voronoi-Creator</b><br />
The Voronoi Diagram creator is located in the voronoi-creator.<br>
Prerequisits:<br>
 *  cmake
 *  libqt4-dev (needs to be tested)
 *  libcgal-dev
 *  libcgal-qt4-dev (needs to be tested, too)
 *  libblas-dev
 *  liblapack-dev
 *  libtbb-dev
 *  libmetis-dev
 *  suitesparse: http://www.cise.ufl.edu/research/sparse/SuiteSparse/

Debian dependencies as one-liner:<br>
`sudo apt-get install cmake libqt4-dev libcgal-dev libcgal-qt4-dev libblas-dev liblapack-dev libtbb-dev libmetis-dev`
<br>
<br>
Do following steps to build voronoi-creator:<br>
`cd voronoi-creator\
mkdir build\
cd build\
qmake ..\
make`<br>
<br>
<b>Import in QT-Creator</b>
 *  Open QT-Creator
 *  Choose File > Open File or Project 
 *  Navigate to voronoi-creation folder
 *  Choose voronoi-creation.pro and confirm Dialog with 'Open'
 *  Specify Build directory/directories (e.g. path-to-repo/voronoi-creation/build)
<br>

<b>Usage</b><br>
Voronoi-Creator expects three arguments:<br>
 1.  Path to image to create voronoi-diagram from
 2.  Output File to write voronoi-diagram to (e.g. voronoi.dat)
 3.  Sites to create (i.e. amount of voronoi-cells)

Example usage:<br>
`./voronoi-creator ~/Pictures/source_irradiance.png voronoi.dat 20000`


### 3. OTM (CGAL)
// Assigned to: Patrick<br>
// Ref: implementation in http://www.geometry.caltech.edu/BlueNoise/<br>
<b>Input</b>
 *  Voronoi Diagram [2]
 *  Singularities [1]
 *  Area Distribution

<b>Output</b>
 *  Power Diagram [3]

### 4. Interpolation-> Natural Neigbors CGAL 
// Assigned to: Cam<br>
<b>Input</b>
 *  Power Diagram [3]

<b>Output</b>
 *  Coordinates (xR) [4]

### 5. Target Surface -> C++ (3D part) 
// Assigned to: Cam<br>
<b>Input</b>
 *  Coordinates (xR) [4]
 *  Target Light Direction (calculated from [4])
 *  Surface Mesh
 *  Incoming Light Direction

<b>Output</b>
 *  3D Mesh (target surface)

