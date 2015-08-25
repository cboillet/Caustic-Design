#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"

#include "glog/logging.h"
#include "SurfaceModel.h"
#include "rendering.h"
#include "utils.h"

class TargetOptimization
{
    vector<glm::vec3> computeNormals;
    vector<glm::vec3> x_sources;
    Model* model;

public:
    TargetOptimization();
    ~TargetOptimization();
    void runOptimization(Model* m, Renderer* renderer);
    void gatherVertexInformation(Vertex * v, uint vertexIndex, vector<int> & neighbors, vector<int> & neighborMap);
    void optimize(Renderer* renderer);
    bool converged();
    void runTest(Renderer* renderer);
};

template<typename T> void cross(T* v1, T* v2, T* result);
template<typename T> void calcFaceNormal(const T* const v1, const T* const v2, const T* const v3, T* result);
template<typename T> T angle(T* v1, T* v2);
template<typename T> void normalize(T* v);
template<typename T> T evaluateInt(const T* const vertex, const T** const neighbors, uint nNeighbors, const vector<int> & neighborMap);
template<typename T> void calcVertexNormal(const T* vertex, T* result, T** faceNormals, const T** neighbors, const vector<int> & neighborMap);
template<typename T> T evaluateReg(const T** const allVertices, const float* L, uint nVertices);

#endif // TARGETOPTIMIZATION_H
