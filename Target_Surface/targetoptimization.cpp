#include "targetoptimization.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class CostFunctorEint{
public:
    CostFunctorEint(int numX, vector<glm::vec3> nt): numX(numX),nT(nt){}

    template <typename T>
    bool operator()(const T* const x1,const T* const x2,const T* const x3, T* e) const{
        //e[0] = T(0);
        for (int i=0; i<numX; i++)
            e[i] = T(i*10.0) - x1[i];
        return true;
    }

private:
    int numX;
    vector<glm::vec3> nT;
};

struct CostFunctorTest {
  template <typename T> bool operator()(const T* const x, T* residual) const {
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
    m.meshes[0].calculateVertexNormals();
    m.setNormals(true);
    for(int i=0; i<m.currentNormals.size(); i++){
        glm::vec3 iN = m.currentNormals[i];
        iN.x = - iN.x;
        m.incidentNormals.push_back(iN);
    }
    m.desiredNormals.reserve(m.meshes[0].selectVerticesMeshFaceNoEdge().size());

    while(!converged(m) &&  numberIteration<1){

        //STEP 1: find light direction
        vector<glm::vec3> Dt = m.computeLightDirectionsScreenSurface();

        //STEP 2: hypothese 1: find desired normals
        m.fresnelMapping();


        //STEP 3: 3D optimization
        optimize(m,Dt);


        //STEP 4: update normals of position
        m.meshes[0].calculateVertexNormals();
        m.setNormals(true);
        m.updateIncidentNormals();

        numberIteration++;
    }
}

void TargetOptimization::optimize(Model& m, vector<glm::vec3> nt){

    int n = NORMALS;

    //no vector<glm::vec3> can be passed to the cost function, but we can run it on each coordinates
    double *x1 = new double[n];
    double *initial_x1 = new double[n];
    double *x2 = new double[n];
    double *initial_x2 = new double[n];
    double *x3 = new double[n];
    double *initial_x3 = new double[n];
    
    

    for (int i=0; i<n; i++)
    {
        x1[i] = m.meshes[0].selectVerticesMeshFaceEdge()[i].Position.x;
        x2[i] = m.meshes[0].selectVerticesMeshFaceEdge()[i].Position.y;
        x3[i] = m.meshes[0].selectVerticesMeshFaceEdge()[i].Position.z;
        initial_x1[i] = x1[i];
        initial_x2[i] = x2[i];
        initial_x3[i] = x3[i];
    }

    Problem problem;

    CostFunction* cost_function_Eint =
        new AutoDiffCostFunction<CostFunctorEint, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint(n, nt));
    //for (int i=0; i<n; i++)

    problem.AddResidualBlock(cost_function_Eint, NULL, x1, x2, x3);

    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

   // std::cout << summary.BriefReport() << std::endl;

//    for (int i=0; i<n; i++)
//        std::cout << "x[" << i << "] : " << initial_x[i] << " -> " << x[i] << std::endl;


    delete[] x1;
    delete[] initial_x1;
    delete[] x2;
    delete[] initial_x2;
    delete[] x3;
    delete[] initial_x3;


}

bool TargetOptimization::converged(Model& m){
    bool converged = true;

    //we check the normals of the face outside the edges

    vector<glm::vec3> vecC = m.currentNormals;
    vector<glm::vec3> vecD = m.desiredNormals;
    std::cout<<"currentNormals size"<<m.currentNormals.size()<<std::endl;
    std::cout<<"desiredNormals size"<<m.desiredNormals.size()<<std::endl;

    for(int i = 0; i<vecC.size(); i++){
        if(glm::distance(vecC[i],vecD[i]) >  CONVERGENCE_LIMIT) converged = false;
    }
    return converged;
}


