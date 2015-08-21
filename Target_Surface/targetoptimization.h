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
    vector<Vertex> x_sources;
    Model* model;

public:
    TargetOptimization();
    ~TargetOptimization();
    void runOptimization(Model* m, Renderer* renderer);
    void optimize(Renderer* renderer);
    bool converged();
};

#endif // TARGETOPTIMIZATION_H
