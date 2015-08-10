#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "SurfaceModel.h"

class TargetOptimization
{
public:
    TargetOptimization();

    void runOptimization(Model& model);
    void runCeresTest();
};

#endif // TARGETOPTIMIZATION_H
