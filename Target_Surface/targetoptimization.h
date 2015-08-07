#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"
#include "glog/logging.h"

class TargetOptimization
{
public:
    TargetOptimization();

    void runCeresTest();
};

#endif // TARGETOPTIMIZATION_H
