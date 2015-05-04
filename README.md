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
 *  Voronoi Diagram [2]

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

