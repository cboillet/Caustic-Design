# Caustic-Design-Project

There are two major parts in this project. 

 1.  Caustic-Design - handles Optimal-Transport and corresponding tasks
 2.  Target Optimization - handles 3D Optimization


## Caustic-Design

Note: This step is optional.. qmake is not fully functional in the current version. Please see the Build-Section for details.

<b>Import in QT-Creator</b>
 *  Open QT-Creator
 *  Choose File > Open File or Project 
 *  Navigate to Caustic_Design folder
 *  Choose Caustic_Design.pro and confirm Dialog with 'Open'
<br>

## Dependencies

Following dependencies are needed for the Caustic_Designer:<br>

 *  cmake
 *  libqt4-dev (needs to be tested)
 *  libcgal-dev
 *  libcgal-qt4-dev (needs to be tested, too)
 *  libblas-dev
 *  liblapack-dev
 *  libtbb-dev
 *  libmetis-dev
 *  libsuitesparse-dev (or manually via: http://faculty.cse.tamu.edu/davis/suitesparse.html)
 *  liblbfgs-dev
 *  libtinyxml-dev

Debian dependencies as one-liner:<br>
`sudo apt-get install cmake libqt4-dev libcgal-dev libcgal-qt4-dev libblas-dev liblapack-dev libtbb-dev libmetis-dev libsuitesparse-dev liblbfgs-dev libtinyxml-dev`

## Build

We suggest using cmake to build the project. To do so, simple:

 1.  Create Build Directory (e.g. `mkdir build-Caustic_Design`)
 2.  Run cmake in the build directory (e.g. `cd build-Caustic_Design && cmake ../Caustic_Design/`)


You may also use qmake instead of cmake if you prefer qmake. But in the current version, suitesparse does not seem to be set correctly when installing it via `apt-get`.


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
algorithm->compute interpolation to load the source image, and points (.dat) and weights (.weights) data from the OTM and run the interpolation

------

## Target Optimization

## Dependencies
Following libraries are needed:

 *  ceres-solver: http://ceres-solver.org/building.html
 *  libassimp-dev

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

