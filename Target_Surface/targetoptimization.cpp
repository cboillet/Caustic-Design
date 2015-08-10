#include "targetoptimization.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

struct CostFunctorTest {
  template <typename T> bool operator()(const T* const x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};

struct CostFunctorEint {
  template <typename T> bool operator()(const T* const x, T* residual, Model& model) const {
    //x position on the surface
    glm::vec3 n = glm::vec3(1,0,0); //incident light ray


    residual[0] = T(10.0) - x[0];
    return true;
  }
};

struct CostFunctorEdir {
  template <typename T> bool operator()(const T* const x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};

struct CostFunctorEflux {
  template <typename T> bool operator()(const T* const x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};

struct CostFunctorEbar {
  template <typename T> bool operator()(const T* x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};


TargetOptimization::TargetOptimization()
{
}

void TargetOptimization::runCeresTest()
{

    double x = 0.5;
    const double initial_x = x;

    Problem problem;

    CostFunction* cost_function =
        new AutoDiffCostFunction<CostFunctorTest, 1, 1>(new CostFunctorTest);
    problem.AddResidualBlock(cost_function, NULL, &x);

    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << "x : " << initial_x
              << " -> " << x << "\n";

}

void TargetOptimization::runOptimization(Model& model){

}
