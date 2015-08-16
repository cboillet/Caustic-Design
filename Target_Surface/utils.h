#ifndef UTILS_H
#define UTILS_H
#include"global.h"
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
#include <math.h>
#include <vector>
#include <string>


using namespace std;

int outTriplet(vector<int> vec, int begin, int end);
bool floatEquals(float val1, float val2);
float fbar(float x, float dth);
glm::vec3 proj(glm::vec3 xs, glm::vec3 di, glm::vec3 pos);
float* matrixProduct(array* L, float* X);
void printMatrix(float* X);
#endif // UTILS_H
