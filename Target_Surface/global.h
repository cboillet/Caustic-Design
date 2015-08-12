#ifndef GLOBAL_H
#define GLOBAL_H
#define GL_GLEXT_PROTOTYPES
#include "GL/glew.h"
#include "GL/glext.h"

#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"
#include "GL/freeglut.h"
#include "GL/glext.h"
#include "glm/glm.hpp"

#include <QGLWidget>
#include <QGLFunctions>

#include <iostream>
#include<numeric>
#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#define MESH_AMOUNT 50 //amount of
#define CAUSTIC_DOMAIN 0.5
#define CONVERGENCE_LIMIT 0.0001
#define NORMALS 1089 //size of vertex in the front face without the edges
#define AIR_REFRACTIV_INDEX 1
#define MATERIAL_REFRACTIV_INDEX 1.49//value for acrylic used in the paper to do test experimently

#endif // GLOBAL_H
