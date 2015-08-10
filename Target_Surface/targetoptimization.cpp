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

struct CostFunctorTest {
  template <typename T> bool operator()(const T* const x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};


class CostFunctorEint{
public:
    CostFunctorEint(){}
    template <typename T>
    bool operator() (const T* x, T* residual) const{
          //x positions on the surface
           for (int i=0; i<4; i++){
           //residual[i] = T(10.0) - x[0];
           }
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


TargetOptimization::TargetOptimization(Model& m)
{
    //computeNormals.reserve(m.meshes[0].selectVerticesMeshFaceEdge().size());
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

}


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
              << " -> " << x[0] << /*", " << x[1] << */ //std::endl;

    //delete[] x;
    //delete[] initial_x;

void TargetOptimization::runOptimization(Model& m){
    while(!converged(m)){
        m.meshes[0].calculateVertexNormals();
        optimize(m);
    }
}

void TargetOptimization::optimize(Model& m){   
    vector<Vertex> x = m.meshes[0].selectVerticesMeshFaceEdge();
    const vector<Vertex> inital_x = x;
    //double y = 0.5;
    //int Eint;
    int para = x.size();
    Problem problem;
    //CostFunction* cost_function = new AutoDiffCostFunction<CostFunctorEint, 1, ceres::DYNAMIC>(new CostFunctorEint());
    //problem.AddResidualBlock(cost_function, NULL, y);
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);


}

bool TargetOptimization::converged(Model& m){
    bool converged = true;
    glm::vec3 norm0;
    glm::vec3 norm1;
    vector<Vertex> vec = m.meshes[0].selectVerticesMeshFaceEdge();
    for(int i = 0; i<m.meshes[0].indices.size(); i++){
        norm0 = vec[i].Normal;
        norm1 = computeNormals[i];
        if(glm::distance(norm0,norm1) >  CONVERGENCE_LIMIT) converged = false;
    }
    return converged;
}

