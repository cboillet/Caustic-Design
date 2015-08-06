#ifndef TARGETOPTIMIZATION_H
#define TARGETOPTIMIZATION_H

#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class TargetOptimization
{
public:
    TargetOptimization();

    void runCeresTest();
};

#endif // TARGETOPTIMIZATION_H
