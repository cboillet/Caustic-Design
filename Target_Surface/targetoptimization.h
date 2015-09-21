#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"

#include "glog/logging.h"
#include "SurfaceModel.h"
#include "utils.h"
#include "rendering.h"
#include "costFunctor.h"

class TargetOptimization
{
    vector<glm::vec3> computeNormals;
    vector<glm::vec3> x_sources;
    Model* model;

public:
    TargetOptimization();
    ~TargetOptimization();
    void runOptimization(Model* m, Renderer* r);
    void gatherVertexInformation(Vertex * v, uint vertexIndex, vector<int> & neighbors, vector<int> & neighborMap, vector<int> & eightNeighbors);
    void addResidualBlocks(ceres::Problem * problem, uint vertexIndex, vector<int> & neighbors, vector<int> & neighborMap, vector<int> & eightNeighbors, double* vertices);

};

#endif // TARGETOPTIMIZATION_H
