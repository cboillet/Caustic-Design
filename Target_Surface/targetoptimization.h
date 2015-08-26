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
    void optimize(Renderer* r);
    void runTest(Renderer* r);
    bool converged();
    void gatherVertexInformation(Vertex * v, uint vertexIndex, vector<int> & neighbors, vector<int> & neighborMap);
    void addResidualBlocks(ceres::Problem * problem, uint vertexIndex, vector<int> & neighbors, vector<int> & neighborMap, double* vertices, Renderer* r);

};

#endif // TARGETOPTIMIZATION_H
