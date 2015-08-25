#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"

#include "glog/logging.h"
#include "SurfaceModel.h"
#include "utils.h"
#include "rendering.h"

class TargetOptimization
{
    vector<glm::vec3> computeNormals;
    vector<Vertex> x_sources;
    Model* model;

public:
    TargetOptimization();
    ~TargetOptimization();
    void runOptimization(Model* m, Renderer* r);
    void optimize(Renderer* r);
    bool converged();
};

#endif // TARGETOPTIMIZATION_H
