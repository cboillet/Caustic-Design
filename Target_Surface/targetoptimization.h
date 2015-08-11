#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "SurfaceModel.h"

class TargetOptimization
{
    vector<glm::vec3> computeNormals;

public:
    TargetOptimization(Model& m);

    void runOptimization(Model& m);
    void runCeresTest();
    void optimize(Model& m, vector<glm::vec3> nt);
    bool converged(Model& m);
};

#endif // TARGETOPTIMIZATION_H
