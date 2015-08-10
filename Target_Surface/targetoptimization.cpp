#include "targetoptimization.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class MyCostFunctor{
public:
    MyCostFunctor(int numX): numX(numX){}

    template <typename T>
    bool operator()(const T* const x, T* e) const{
        for (int i=0; i<numX; i++)
        {
            e[0] += T(10.0) - x[i];
        }
        return true;
    }

private:
    int numX;
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

/*
 double *x = new double[1];
    x[0] = 3;
    //x[1] = 2;
    double *initial_x = new double[1];
    initial_x[0] = 3;

    Problem problem;

    CostFunction* cost_function =
        new AutoDiffCostFunction<MyCostFunctor, 1, ceres::DYNAMIC>(new MyCostFunctor(1));
    problem.AddResidualBlock(cost_function, NULL, x);

    // Run the solver!
    Solver::Options options;
    //options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << "x[0] : " << initial_x[0] //<< ", x[1]:" << initial_x[1]
              << " -> " << x[0] << /*", " << x[1] << */std::endl;

    //delete[] x;
    //delete[] initial_x;
*/
}

void TargetOptimization::runOptimization(Model& model){

}

