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
    CostFunctorEint(Model* m, float w): model(m), weight(w){
        vertices = m->meshes[0].faceVertices;
    }

    bool operator()(const double* x1, const double* x2, const double* x3, double* e) const{
        //load the positions
        for (int i=0; i<NORMALS; i++){
            vertices[i]->Position.x=x1[i];
            vertices[i]->Position.y=x2[i];
            vertices[i]->Position.z=x3[i];
        }
        model->meshes[0].calculateVertexNormals();
        //model->computeLightDirectionsScreenSurface();
        //model->fresnelMapping();

        for (int i=0; i<NORMALS; i++)
                //if(!model->meshes[0].isEdge(vertices[i])){
                e[i] = weight*(glm::length(vertices[i]->Normal-model->desiredNormals[i]));
           // }
           // else e[i]=0;
        return true;
    }

private:
    Model* model;
    vector<Vertex*> vertices;
    float weight;
};


class CostFunctorEdir{
public:
    CostFunctorEdir(vector<Vertex> s, float w): sources(s), weight(w){
    }

    template <typename T> bool operator()(const T* const x1, const T* const x2, const T* const x3, T* e) const{
        for (int i=0; i<NORMALS; i++){
            // x - proj(..) will only return the difference in y- and z-coordinates. we may ignore the rest
            T y = x2[i] - T(sources[i].Position.y);
            T z = x3[i] - T(sources[i].Position.z);
            e[i] = T(weight) * (ceres::abs(y)+ceres::abs(z));
        }
        return true;
    }

private:
    vector<Vertex> sources;
    float weight;
};

class CostFunctorTest {
  public:
    CostFunctorTest(Model* m):model(m){
                vector<Vertex*> vertices = m->meshes[0].faceVertices;
    }
    template <typename T> bool operator()(const T* const x1, T* e) const {
        vertices[0]->Position.x=x1[0];
                for (int i=0; i<NORMALS; i++)
                    e[i]=T(10.0) - x1[0];
    return true;
  }
private:
    Model* model;
    vector<Vertex*> vertices;
};


struct CostFunctorEflux {
  template <typename T> bool operator()(const T* const x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};


class CostFunctorEbar {
public:
    CostFunctorEbar (Model* m, float w): model(m), weight(w)
    {
        for (int i=0; i<NORMALS; i++)
        {
            receiverPositions.push_back(&model->receiverLightPositions[i]);
        }

        surfaceVertices = model->meshes[0].faceVertices;
    }


    template <typename T> bool operator()(const T* const x1, T* e) const{

        T dth = T(EBAR_DETH);

        for (int i=0; i<NORMALS; i++){
            // nr * (x - xr) means dot product with a normal (1,0,0). which means that y- and z-values are not included anyways.
            // so we only substract the x-positions (and dot product would then only multiply with 1)
            T x = x1[i] - T(receiverPositions[i]->x);
            e[i] = T(weight) * ceres::max(T(0), - ceres::log( (T(1)-x) + dth) );
        }

        return true;
    }



private:
    Model* model;
    float weight;
    vector<glm::vec3*> receiverPositions;
    vector<Vertex*> surfaceVertices;
};


class CostFunctorEreg{
public:
    CostFunctorEreg(Model* m, float w): model(m), weight(w){
        vertices = m->meshes[0].faceVertices;
        for (int i=0; i<NORMALS; i++){
                for (int j=0; j<NORMALS; j++){
                    L[i][j]=0;
                }

                vector<int> neighbors = model->meshes[0].getNeighborsIndex(vertices[i]);
                //vector<int> neighbors = model->meshes[0].getClosestNeighbors(vertices[i]);
                //std::cout<<"neighbors.size"<<neighbors.size()<<std::endl;
                L[i][i]=-neighbors.size();
                //std::cout<<"L[i][i]"<<L[i][i]<<std::endl;
                for(int j=0; j<neighbors.size(); j++){
                    if(j!=i)
                        L[i][neighbors[j]]=1;
                }
                neighbors.clear();
        }
    }

    ~CostFunctorEreg(){
        model = NULL;
        vertices.clear();
        delete[] L;
    }

    bool operator()(const double* x1, const double* x2, const double* x3, double* e) const{
        double result;
        for(int i=0; i<NORMALS; i++){
            vertices[i]->Position.x = x1[i];
            vertices[i]->Position.y = x2[i];
            vertices[i]->Position.z = x3[i];
        }

        //printMatrix(L);
        for(int i=0; i<NORMALS; i++){
            result=0;
            for (int j=0; j<NORMALS; j++){
                result +=
                        (L[i][j]*(vertices[j]->Position.x)) +
                        (L[i][j]*(vertices[j]->Position.y)) +
                        (L[i][j]*(vertices[j]->Position.z));

                //result += weight*(pow(L[i][j]*(vertices[j]->Position.x),2)+pow(L[i][j]*(vertices[j]->Position.y),2)+pow(L[i][j]*vertices[j]->Position.z,2));
            }
            e[i] = weight * fabs(result);
        }
        return true;
    }

private:
    Model* model;
    vector<Vertex*> vertices;
    float weight;
    array* L = new array[NORMALS];
};


class CostFunctorEreg2{
public:
    CostFunctorEreg2(Model* m, float w): model(m), weight(w){
        for (int i=0; i<NORMALS; i++){
                for (int j=0; j<NORMALS; j++){
                    L[i][j]=0;
                }

                vector<int> neighbors = model->meshes[0].getNeighborsIndex(model->meshes[0].faceVertices[i]);

                L[i][i]= - int(neighbors.size());
                for(int j=0; j<neighbors.size(); j++){
                    if(neighbors[j] != i)
                    {
                        L[i][neighbors[j]] = 1;
                    }
                }
                neighbors.clear();
        }
    }

    template <typename T> bool operator()(const T* const x1,const T* const x2,const T* const x3, T* e) const {
        T res1, res2, res3;

        for(int i=0; i<NORMALS; i++){
            res1 = res2 = res3 = T(0);
            for (int j=0; j<NORMALS; j++){
                res1 += T(L[i][j]) * x1[j];
                res2 += T(L[i][j]) * x2[j];
                res3 += T(L[i][j]) * x3[j];
            }

            e[i] = T(weight) * (ceres::abs(res1) + ceres::abs(res2) + ceres::abs(res3));
        }
        return true;
    }

private:
    Model* model;
    float weight;
    array* L = new array[NORMALS];
};

TargetOptimization::TargetOptimization()
{
    //computeNormals.reserve(m.meshes[0].selectVerticesMeshFaceEdge().size());
}

TargetOptimization::~TargetOptimization(){
    if (model != NULL) delete[] model;
}




void TargetOptimization::runOptimization(Model* m){
    int numberIteration = 0;
    model=m;
    for(int i=0; i<m->meshes[0].faceVertices.size(); i++){
        Vertex v = *(m->meshes[0].faceVertices[i]);
        x_sources.push_back(v);
    }
    //load the values of the calculated vertex in vector<glm::vec3> currentNormals we take the edges
    model->meshes[0].calculateVertexNormals();
    model->setNormals(true); //set current Normals

    //for testing the convergence we compute the desired normals one time first
    model->fresnelMapping();


    while(!converged() &&  numberIteration<1){

        //STEP 1: Compute light direction
        model->computeLightDirectionsScreenSurface();
        //STEP 2: update normals of position
//        model->meshes[0].calculateVertexNormals();
//        model->setNormals(true);
        model->fresnelMapping();

        //STEP 1: 3D optimization
        optimize();

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
        x2[i] = x[i]->Position.y;
        x3[i] = x[i]->Position.z;
        initial_x1[i] = x1[i];
        initial_x2[i] = x2[i];
        initial_x3[i] = x3[i];
    }

    Problem problem;

    /*1. Test passing vector<Vertex*>*/
   CostFunction* cost_function_Eint =
        new NumericDiffCostFunction<CostFunctorEint, ceres::CENTRAL ,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint(model, EINT_WEIGHT));
   CostFunction* cost_function_Ebar =
           new AutoDiffCostFunction<CostFunctorEbar,NORMALS, NORMALS>(new CostFunctorEbar(model,EBAR_WEIGHT));
      CostFunction* cost_function_Edir =
           new AutoDiffCostFunction<CostFunctorEdir, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEdir(x_sources, float(EDIR_WEIGHT)));
      //CostFunction* cost_function_Ereg =
      //     new NumericDiffCostFunction<CostFunctorEreg, ceres::CENTRAL ,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEreg(model, EREG_WEIGHT));
      CostFunction* cost_function_Ereg2 =
           new AutoDiffCostFunction<CostFunctorEreg2, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEreg2(model, EREG_WEIGHT));
//   CostFunction* cost_function_Etest =
//        new AutoDiffCostFunction<CostFunctorTest, NORMALS, NORMALS>(new CostFunctorTest(model));


   problem.AddResidualBlock(cost_function_Eint, NULL, x1, x2, x3);
   problem.AddResidualBlock(cost_function_Ebar, NULL, x1);
   problem.AddResidualBlock(cost_function_Edir, NULL, x1, x2, x3);
   problem.AddResidualBlock(cost_function_Ereg2, NULL, x1, x2, x3);
//   problem.AddResidualBlock(cost_function_Ereg2, NULL, x1, x2, x3);
//  problem.AddResidualBlock(cost_function_Etest, NULL, x1);

   /*2. passing template method*/
//    CostFunction* cost_function_Eint2 =
//         new AutoDiffCostFunction<CostFunctorEint2, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint2(model));
//    problem.AddResidualBlock(cost_function_Eint2, NULL, x1, x2, x3);

    /*1. passing template method*/
    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::ITERATIVE_SCHUR; //large bundle adjustment problems
    options.max_num_iterations = 150;
    options.dense_linear_algebra_library_type = ceres::LAPACK;
    options.visibility_clustering_type = ceres::SINGLE_LINKAGE;
    //options.preconditioner_type = ceres::CLUSTER_TRIDIAGONAL; // fast preconditioner
    string error;
    if(!options.IsValid(&error))
    {
        std::cout << "Options not valid: " << error << std::endl;
    }

    Solver::Summary summary;
    Solve(options, &problem, &summary);

   // std::cout << summary.BriefReport() << std::endl;
    std::cout << summary.FullReport() << std::endl;


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

    vector<Vertex*>  vecC = model->meshes[0].faceVertices;
    vector<glm::vec3> vecD = model->desiredNormals;
    std::cout<<"faceVertices size"<<vecC.size()<<std::endl;
    std::cout<<"desiredNormals size"<<model->desiredNormals.size()<<std::endl;

    for(int i = 0; i<vecC.size(); i++){
        if(glm::distance(vecC[i]->Normal,vecD[i]) >  CONVERGENCE_LIMIT) converged = false;
    }
    return converged;
}


