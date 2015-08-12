#include "targetoptimization.h"

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class CostFunctorEint{
public:
    CostFunctorEint(Model* m): model(m){}

    bool operator()(const double* const x1, const double* const x2, const double* const x3, double* e) const{
        vector<Vertex*> surfaceVertices = model->meshes[0].faceVertices;
        //load the positions
        for (int i=0; i<NORMALS; i++){
//            surfaceVertices[i]->Position.x=(float)x1;
//            surfaceVertices[i]->Position.x=x2;
//            surfaceVertices[i]->Position.x=x3;
        }
        model->meshes[0].calculateVertexNormals();
        //e[0] = T(0);
        for (int i=0; i<NORMALS; i++)
            e[i] = glm::length(surfaceVertices[i]->Normal-model->desiredNormals[i]);
        return true;
    }

private:
    Model* model;
};

class CostFunctorEint2{
public:
    CostFunctorEint2(Model* m): model(m){}

    template <typename T>
    bool operator()(const T* const x1,const T* const x2,const T* const x3, T* e) const{
        e[0] = T(10.0) - x1[0];
        return true;
    }

private:
    Model* model;
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


class CostFunctorEbar {
public:
    CostFunctorEbar (Model* m): model(m){}

    template <typename T>
    bool operator()(const T* const x1, T* e) const{
        vector<Vertex*> surfaceVertices = model->meshes[0].faceVertices;
        float dth= model->getFocalLength() + model->meshes[0].getMaxX();
        //load the positions
        for (int i=0; i<NORMALS; i++){
            surfaceVertices[i]->Position.x=x1;
            //normal to the receiver plane
            glm::vec3 nr;
            nr.x = 1;
            nr = 0;
            nr.z = 0;
            e[i] = fbar(glm::dot(nr,(surfaceVertices[i]->Position-model->receiverLightPositions[i])),dth);
        }

        return true;
    }

private:
    Model* model;
};



TargetOptimization::TargetOptimization()
{
    //computeNormals.reserve(m.meshes[0].selectVerticesMeshFaceEdge().size());
}

TargetOptimization::~TargetOptimization(){
    if (model != NULL) delete[] model;
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


void TargetOptimization::runOptimization(Model* m){
    int numberIteration = 0;
    model=m;
    //load the values of the calculated vertex in vector<glm::vec3> currentNormals we take the edges
    model->meshes[0].calculateVertexNormals();
    model->setNormals(true);
    for(int i=0; i<model->currentNormals.size(); i++){
        glm::vec3 iN = model->currentNormals[i];
        iN.x = - iN.x;
        model->incidentNormals.push_back(iN);
    }

    //STEP 2: hypothese 1: find desired normals
    model->fresnelMapping();


    while(!converged() &&  numberIteration<1){
        //STEP 1: 3D optimization
        optimize();


        //STEP 2: update normals of position
        model->meshes[0].calculateVertexNormals();
        model->setNormals(true);
        model->computeLightDirectionsScreenSurface();


        //STEP 4: hypothese 1: find desired normals
        model->fresnelMapping();

        numberIteration++;
    }
}

void TargetOptimization::optimize(){

    int n = NORMALS;

    //no vector<glm::vec3> can be passed to the cost function, but we can run it on each coordinates
    double *x1 = new double[n];
    double *initial_x1 = new double[n];
    double *x2 = new double[n];
    double *initial_x2 = new double[n];
    double *x3 = new double[n];
    double *initial_x3 = new double[n];


    vector<Vertex*> x=model->meshes[0].faceVertices;

    for (int i=0; i<n; i++)
    {
        x1[i] = x[i]->Position.x;
        x1[i] = x[i]->Position.x;
        x1[i] = x[i]->Position.x;
        initial_x1[i] = x1[i];
        initial_x2[i] = x2[i];
        initial_x3[i] = x3[i];
    }

    Problem problem;

    /*1. Test passing vector<Vertex*>*/
//   CostFunction* cost_function_Eint =
//        new NumericDiffCostFunction<CostFunctorEint, ceres::CENTRAL ,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint(model));
//    problem.AddResidualBlock(cost_function_Eint, NULL, x1, x2, x3);

   /*2. passing template method*/
    CostFunction* cost_function_Eint2 =
         new AutoDiffCostFunction<CostFunctorEint2, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint2(model));
    problem.AddResidualBlock(cost_function_Eint2, NULL, x1, x2, x3);

    /*1. passing template method*/
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

bool TargetOptimization::converged(){
    bool converged = true;

    //we check the normals of the face outside the edges

    vector<glm::vec3> vecC = model->currentNormals;
    vector<glm::vec3> vecD = model->desiredNormals;
    std::cout<<"currentNormals size"<<model->currentNormals.size()<<std::endl;
    std::cout<<"desiredNormals size"<<model->desiredNormals.size()<<std::endl;

    for(int i = 0; i<vecC.size(); i++){
        if(glm::distance(vecC[i],vecD[i]) >  CONVERGENCE_LIMIT) converged = false;
    }
    return converged;
}


