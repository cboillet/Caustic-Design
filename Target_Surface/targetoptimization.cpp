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
    CostFunctorEdir(Model* m, vector<Vertex> s, float w): model(m), sources(s), weight(w){
        vertices = m->meshes[0].faceVertices;
    }

    bool operator()(const double* x1, const double* x2, const double* x3, double* e) const{
        //load the positions
        for (int i=0; i<NORMALS; i++){
            vertices[i]->Position.x=x1[i];
            vertices[i]->Position.y=x2[i];
            vertices[i]->Position.z=x3[i];
        }

        for (int i=0; i<NORMALS; i++)
                {
                e[i] = weight*(glm::length(vertices[i]->Position - proj(sources[i].Position, glm::vec3(1,0,0), vertices[i]->Position)));
            }
        return true;
    }

private:
    Model* model;
    vector<Vertex*> vertices;
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
    CostFunctorEbar (Model* m, float w): model(m), weight(w){}
    bool operator()(const double* x1,double* e) const{
        vector<Vertex*> surfaceVertices = model->meshes[0].faceVertices;
        float dth= EBAR_DETH; // model->getFocalLength() + model->meshes[0].getMaxX();
        int j=0;
        //load the positions
        for (int i=0; i<NORMALS; i++){
            surfaceVertices[i]->Position.x=x1[i];
            //normal to the receiver plane
            glm::vec3 nr;
            nr.x = 1;
            nr.y = 0;
            nr.z = 0;
            //if(!model->meshes[0].isEdge(surfaceVerticesEdge[i])){
                e[i] = weight*(fbar(glm::dot(nr,(surfaceVertices[i]->Position-model->receiverLightPositions[i])),dth));
                j++;
            //}
            //else e[i] = 0;
        }

        return true;
    }

private:
    Model* model;
    float weight;
};


class CostFunctorEreg{
public:
    CostFunctorEreg(Model* m, float w): model(m), weight(w){
        vertices = m->meshes[0].faceVertices;
    }

    ~CostFunctorEreg(){
        model = NULL;
        vertices.clear();
        delete[] L;
    }

    bool operator()(const double* x1, const double* x2, const double* x3, double* e) const{
        for(int i=0; i<NORMALS; i++){
            vertices[i]->Position.x = x1[i];
            vertices[i]->Position.y = x2[i];
            vertices[i]->Position.z = x3[i];
        }
        for (int i=0; i<NORMALS; i++){
                for (int j=0; j<NORMALS; j++){
                    L[i][j]=0;
                }

                vector<int> neighbors = model->meshes[0].getNeighborsIndex(vertices[i]);
                L[i][i]=neighbors.size();

                for(int j=0; j<neighbors.size(); j++){
                    L[i][neighbors[j]]=1;
                }
                neighbors.clear();
        }
        for(int i=0; i<NORMALS; i++){
            for (int j=0; j<NORMALS; j++){
                e[i] += weight*(pow(L[i][j]*(vertices[j]->Position.x),2)+pow(L[i][j]*(vertices[j]->Position.y),2)+pow(L[i][j]*vertices[j]->Position.z,2));
            }

        }
        return true;
    }

private:
    Model* model;
    vector<Vertex*> vertices;
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
        new NumericDiffCostFunction<CostFunctorEbar, ceres::CENTRAL ,NORMALS, NORMALS>(new CostFunctorEbar(model,EBAR_WEIGHT));
   CostFunction* cost_function_Edir =
        new NumericDiffCostFunction<CostFunctorEdir, ceres::CENTRAL ,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEdir(model, x_sources, EDIR_WEIGHT));
   CostFunction* cost_function_Ereg =
        new NumericDiffCostFunction<CostFunctorEreg, ceres::CENTRAL ,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEreg(model, EREG_WEIGHT));

//   CostFunction* cost_function_Etest =
//        new AutoDiffCostFunction<CostFunctorTest, NORMALS, NORMALS>(new CostFunctorTest(model));


   problem.AddResidualBlock(cost_function_Eint, NULL, x1, x2, x3);
   problem.AddResidualBlock(cost_function_Ebar, NULL, x1);
   problem.AddResidualBlock(cost_function_Edir, NULL, x1, x2, x3);
   problem.AddResidualBlock(cost_function_Ereg, NULL, x1, x2, x3);
//   problem.AddResidualBlock(cost_function_Etest, NULL, x1);

   /*2. passing template method*/
//    CostFunction* cost_function_Eint2 =
//         new AutoDiffCostFunction<CostFunctorEint2, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint2(model));
//    problem.AddResidualBlock(cost_function_Eint2, NULL, x1, x2, x3);

    /*1. passing template method*/
    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;

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


