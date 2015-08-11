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
        //e[0] = T(0);
        for (int i=0; i<numX; i++)
            e[i] = T(i*10.0) - x[i];
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


void TargetOptimization::runOptimization(Model& m){
    int numberIteration = 0;

    //load the values of the calculated vertex in vector<glm::vec3> currentNormals we take the edges
    m.setNormals(true);
    m.meshes[0].calculateVertexNormals();

    while(!converged(m) &&  numberIteration<1){
        optimize(m);
        m.meshes[0].calculateVertexNormals();
        m.setNormals(true);
        numberIteration++;
    }
}

void TargetOptimization::optimize(Model& m){
//    vector<Vertex> x = m.meshes[0].selectVerticesMeshFaceEdge();
//    const vector<Vertex> inital_x = x;
//    double x = 0.5;
//    const double initial_x = x;
//    double x[] = {1.0, 2.0, 3.0};
//    Problem problem;
//    CostFunction* cost_function =
//        new AutoDiffCostFunction<CostFunctorEint, 1, 1>(new CostFunctorEint);
//    problem.AddResidualBlock(cost_function, NULL, x);
//    Solver::Options options;
//    options.minimizer_progress_to_stdout = true;
//    Solver::Summary summary;
//    Solve(options, &problem, &summary);

    int n = NORMALS;
    double *x = new double[n];
    double *initial_x = new double[n];

    for (int i=0; i<n; i++)
    {
        x[i] = i+1;
        initial_x[i] = i+1;
    }
    Problem problem;

    CostFunction* cost_function =
        new AutoDiffCostFunction<MyCostFunctor, ceres::DYNAMIC, NORMALS>(new MyCostFunctor(n), n);
    //for (int i=0; i<n; i++)
        problem.AddResidualBlock(cost_function, NULL, x);

    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

   // std::cout << summary.BriefReport() << std::endl;

    for (int i=0; i<n; i++)
        std::cout << "x[" << i << "] : " << initial_x[i] << " -> " << x[i] << std::endl;

    delete[] x;
    delete[] initial_x;


}

bool TargetOptimization::converged(Model& m){
    bool converged = true;
    vector<glm::vec3> vecC = m.currentNormals;
    vector<glm::vec3> vecD = m.desiredNormals;
    std::cout<<"currentNormals size"<<m.currentNormals.size()<<std::endl;
    std::cout<<"desiredNormals size"<<m.desiredNormals.size()<<std::endl;

    for(int i = 0; i<vecC.size(); i++){
        if(glm::distance(vecC[i],vecD[i]) >  CONVERGENCE_LIMIT) converged = false;
    }
    return converged;
}

