# Caustic-Design

<b>Import in QT-Creator</b>
 *  Open QT-Creator
 *  Choose File > Open File or Project 
 *  Navigate to Caustic_Design folder
 *  Choose Caustic_Design.pro and confirm Dialog with 'Open'
<br>

## Dependencies

Following dependencies are needed for this branch:

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
 *  libtinyxml-dev

Debian dependencies as one-liner:<br>
`sudo apt-get install cmake libqt4-dev libcgal-dev libcgal-qt4-dev libblas-dev liblapack-dev libtbb-dev libmetis-dev libtinyxml-dev`

## Usage

Following sections describe the usage of the different steps for the caustic-design project.

### Optimal Transport

To run the optimal transport, a source-image as well as a target-image needed. The rest of the steps are done automatically.<br>

----------

<b>Loading source- and target-image</b><br>
 *  Via code:  In window.cpp uncomment the lines `//open(QString("/home/p/Pictures/einstein_small.png"), false);` and `//open(QString("/home/p/Pictures/white_small.png"), true);` and replace the strings with the path to the target image and source image, respectively.
 *  Via UI: `File > Load Image` to load target image and `File > Load Source Image` to load source image.

<b>NOTE</b>: Source- and Target-Image need to be of same ratio.

---------

<b>Configuration</b><br>
All relevant values (amount of sites, multi-scale levels) can be configured in the config.h file.

---------

<b>Running</b>
To finally run the Optimal Transport: 
 *  `Algorithm > Compute Optimal Transport`
 *  When running Gradient Descent, ensure that `LEVEL_MAX` is set to `1`


------

### Interpolation-> Natural Neigbors CGAL 
// Assigned to: Cam<br>
To run the interpolation after computing the Optimal transport, load on the left scene (m_scene) the source irradiance, compute the interpolation, results appear on the left scene.

------


### Target Surface -> C++ (3D part) 
// Assigned to: Cam<br>
<b>Input</b>
 *  Coordinates (xR) [4]
 *  Target Light Direction (calculated from [4])
 *  Surface Mesh
 *  Incoming Light Direction

<b>Output</b>
 *  3D Mesh (target surface) 

<b>Computing the surface optomization</b><br />
The code is located in target-surface-optimization folder<br>
Prerequisits:<br>
 *  openGL
 *  GLM header librairy
 *  glew
 *  glfw
 * SOIL
 * assimp

Debian dependencies as one-liner:<br>
`sudo apt-get install libglew-dev libsoil-dev`
<br>

