#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"

#include "glog/logging.h"
#include "SurfaceModel.h"
#include "rendering.h"
#include "utils.h"
#include "cost_functors.h"

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

#endif // TARGETOPTIMIZATION_H
