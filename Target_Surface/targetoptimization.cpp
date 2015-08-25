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

    bool operator()(const double* x1/*, const double* x2, const double* x3*/, double* e) const{
        //load the positions
        for (int i=0; i<NORMALS; i++){
            vertices[i]->Position.x=x1[i];
            //vertices[i]->Position.y=x2[i];
            //vertices[i]->Position.z=x3[i];
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

class CostFunctorEint2{
public:
    CostFunctorEint2(Model* m, float w): model(m), weight(w){
        allVertices = m->meshes[0].allVertices;
        vertices = m->meshes[0].faceVertices;
        edgeVertices = m->meshes[0].faceVerticesEdge;
        edgeIndices = m->meshes[0].edgeIndices;
        edgeNoEdgeMapping = m->meshes[0].edgeToNoEdgeMapping; //createEdgeToNoEdgeMapping();
        //noEdgeEdgeMapping = m->meshes[0].noEdgeToEdgeMapping; //createNoEdgeToEdgeMapping();
        edgeAdjacentFaces = m->meshes[0].edgeAdjacentFaces;
        //compareVectors(edgeNoEdgeMapping , m->meshes[0].edgeToNoEdgeMapping);
    }

    template <typename T> bool operator()(const T* const x1, const T* const x2, const T* const x3, T* e) const{

        // --------- face normals -----------
        T* v1 = new T[3];
        T* v2 = new T[3];

        T** faceNormals  = new T*[edgeIndices.size()];
        for (uint i=0; i<edgeIndices.size(); i++)
        {
            calcFaceNormal(faceNormals, v1, v2, x1, x2, x3, i);
        }

        // ----------- vertex normals ------------

        T** vertexNormals = new T*[edgeVertices.size()];
        for (uint i=0; i<edgeVertices.size(); i++)
        {
            calcVertexNormal(vertexNormals, faceNormals, x1, x2, x3, i);
        }

        // ----------- evaluation -------------


        for (uint i=0; i<edgeVertices.size(); i++)
        {
            int index = edgeNoEdgeMapping[i];

            if(index != -1)
            {
                T x = T(model->desiredNormals[index].x) - vertexNormals[i][0];
                T y = T(model->desiredNormals[index].y) - vertexNormals[i][1];
                T z = T(model->desiredNormals[index].z) - vertexNormals[i][2];
                T res = ceres::sqrt(x*x + y*y + z*z);
                e[index] = res;

            }
        }



        // ------------ clear ------------

        delete[] v1;
        delete[] v2;

        for(uint i=0; i<edgeIndices.size(); i++)
        {
            delete[] faceNormals[i];

            if(i < vertices.size())
                delete[] vertexNormals[i];
        }


        delete[] faceNormals;
        delete[] vertexNormals;

        // ------------ done --------------

        std::cerr << "done " << std::endl;

        return true;
    }

    template<typename T> void calcFaceNormal(T** faceNormals, T* v1, T* v2, const T* const x1, const T* const x2, const T* const x3, uint i) const
    {
        int index = edgeIndices[i][0];
        T vertex1[3];
        assign2(&vertex1[0], index, x1, x2, x3);

        index = edgeIndices[i][1];
        T vertex2[3];
        assign2(&vertex2[0], index, x1, x2, x3);

        index = edgeIndices[i][2];
        T vertex3[3];
        assign2(&vertex3[0], index, x1, x2, x3);

        v1[0] = vertex2[0] - vertex1[0];
        v1[1] = vertex2[1] - vertex1[1];
        v1[2] = vertex2[2] - vertex1[2];

        v2[0] = vertex3[0] - vertex1[0];
        v2[1] = vertex3[1] - vertex1[1];
        v2[2] = vertex3[2] - vertex1[2];

        T* result = new T[3];
        cross(v1, v2, result);

         faceNormals[i] = result;
    }

    template<typename T> void calcVertexNormal(T** vertexNormals, T** faceNormals, const T* const x1, const T* const x2, const T* const x3, uint index) const
    {
        T* result = new T[3];
        result[0] = result[1] = result[2] = T(0);

        T edge1Res[3];
        T edge2Res[3];
        T edge1[3];
        T edge2[3];
        T edge1Sub[3];
        T edge2Sub[3];

        for(uint i=0; i<edgeAdjacentFaces[index].size(); i++)
        {
            // find out which vertex of the current face is the vertex we are currently looking at
            // aF[vertexIndex] is a list of faces (aka a list of indices of the indices-vector)
            // so indices[aF[vertexIndex][j]] is a glm::vec3 that contains one face
            // and the current vertex is vertex[vertexIndex]
            int thisVertexIndex = -1;
            for(int vIndex=0; vIndex < 3; vIndex++)
            {

                int indx = edgeIndices[edgeAdjacentFaces[index][i]][vIndex];
                T pos1[3];
                assign2(&pos1[0], indx, x1, x2, x3);


                T pos2[3];
                assign2(&pos2[0], index, x1, x2, x3);

                if(ceres::sqrt(
                               // x²
                               (pos1[0]-pos2[0])*(pos1[0]-pos2[0]) +
                               // y²
                               (pos1[1]-pos2[1])*(pos1[1]-pos2[1]) +
                               // z²
                               (pos1[2]-pos2[2])*(pos1[2]-pos2[2]))
                                < T(0.0001))
                {
                    thisVertexIndex = vIndex;
                    break;
                }

            }

            // we got index of our current vertex within the face, now get others
            int other1 = (thisVertexIndex+1) % 3;
            int other2 = (thisVertexIndex+2) % 3;

            // create the vectors the represent the edges from current vertex to the other 2
            int indx = edgeIndices[edgeAdjacentFaces[index][i]][thisVertexIndex];
            assign2(&edge1[0], indx, x1, x2, x3);

            indx = edgeIndices[edgeAdjacentFaces[index][i]][other1];
            assign2(&edge1Sub[0], indx, x1, x2, x3);


            indx = edgeIndices[edgeAdjacentFaces[index][i]][thisVertexIndex];
            assign2(&edge2[0], indx, x1, x2, x3);

            indx = edgeIndices[edgeAdjacentFaces[index][i]][other2];
            assign2(&edge2Sub[0], indx, x1, x2, x3);

            for(uint j=0; j<3; j++)
            {
                edge1Res[j] = edge1[j] - edge1Sub[j];
                edge2Res[j] = edge2[j] - edge2Sub[j];
            }

            normalize(edge1);
            normalize(edge2);

            T incidentAngle = angle(edge1Res, edge2Res);

            /*// get angle between the edges
            float incidentAngle = abs(glm::angle(glm::normalize(edge1), glm::normalize(edge2)));
            if(incidentAngle > 180)
               incidentAngle = 360 - incidentAngle;*/

            // use that angle as weighting
            //vertexNormal += (faceNormals[adjacentFaces[vertexIndex][j]] * incidentAngle);
            for (uint j=0; j<3; j++)
            {
                result[j] += faceNormals[edgeAdjacentFaces[index][i]][j] * incidentAngle;
            }

        }


        vertexNormals[index] = result;

        normalize(vertexNormals[index]);
    }

    template<typename T> T angle(T* v1, T* v2) const
    {
        return ceres::acos(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    }

    template<typename T> void normalize(T* v) const
    {
        T sum = T(0);
        for (uint i=0; i<3; i++)
        {
            sum += v[i];
        }

        for (uint i=0; i<3; i++)
        {
            v[i] = v[i]/sum;
        }
    }


    template<typename T> void cross(T* v1, T* v2, T* result) const
    {
        result[0] = T(v1[1]*v2[2] - v1[2]*v2[1]);
        result[1] = T(v1[2]*v2[0] - v1[0]*v2[2]);
        result[2] = T(v1[0]*v2[1] - v1[1]*v2[0]);
    }

    template<typename T> void assign( T (&to)[3], int index, const T* const x1, const T* const x2, const T* const x3) const
    {

        int mappedIndex = edgeNoEdgeMapping[index];

        if(mappedIndex == -1)
        {
            to[0] = T(edgeVertices[index]->Position.x);
            to[1] = T(edgeVertices[index]->Position.y);
            to[2] = T(edgeVertices[index]->Position.z);
        }else
        {
            to[0] = x1[mappedIndex];
            to[1] = x2[mappedIndex];
            to[2] = x3[mappedIndex];
        }
    }

    template<typename T> void assign2( T* to, int index, const T* const x1, const T* const x2, const T* const x3) const
    {

        int mappedIndex = edgeNoEdgeMapping[index];

        if(mappedIndex == -1)
        {
            to[0] = T(edgeVertices[index]->Position.x);
            to[1] = T(edgeVertices[index]->Position.y);
            to[2] = T(edgeVertices[index]->Position.z);
        }else
        {
            to[0] = x1[mappedIndex];
            to[1] = x2[mappedIndex];
            to[2] = x3[mappedIndex];
        }
    }

private:
    Model* model;
    vector<Vertex*> allVertices;
    vector<Vertex*> vertices;
    vector<Vertex*> edgeVertices;
    vector<int> edgeNoEdgeMapping;
    vector<int> noEdgeEdgeMapping;
    vector<glm::uvec3> edgeIndices;
    vector<vector<uint> > edgeAdjacentFaces;
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
    CostFunctorEreg2(Model* m, Renderer* renderer, float w): model(m), weight(w){
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
        //printMatrix(L);
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

            e[i] = T(weight) * res1 + ceres::abs(res2) + ceres::abs(res3);
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




void TargetOptimization::runOptimization(Model* m, Renderer* renderer){
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
        optimize(renderer);

        numberIteration++;
    }
}

void TargetOptimization::optimize(Renderer* renderer){

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
        new AutoDiffCostFunction<CostFunctorEint2,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEint2(model, EINT_WEIGHT));

           CostFunction* cost_function_Ebar =
           new AutoDiffCostFunction<CostFunctorEbar,NORMALS, NORMALS>(new CostFunctorEbar(model,EBAR_WEIGHT));
      //CostFunction* cost_function_Edir =
      //     new AutoDiffCostFunction<CostFunctorEdir, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEdir(x_sources, float(EDIR_WEIGHT)));
      //CostFunction* cost_function_Ereg =
      //     new NumericDiffCostFunction<CostFunctorEreg, ceres::CENTRAL ,NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEreg(model, EREG_WEIGHT));
      CostFunction* cost_function_Ereg2 =
           new AutoDiffCostFunction<CostFunctorEreg2, NORMALS, NORMALS, NORMALS, NORMALS>(new CostFunctorEreg2(model, renderer, EREG_WEIGHT));
//   CostFunction* cost_function_Etest =
//        new AutoDiffCostFunction<CostFunctorTest, NORMALS, NORMALS>(new CostFunctorTest(model));


   problem.AddResidualBlock(cost_function_Eint, NULL, x1, x2, x3);
   problem.AddResidualBlock(cost_function_Ebar, NULL, x1);
   //problem.AddResidualBlock(cost_function_Edir, NULL, x1, x2, x3);
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
    options.max_num_iterations = 200;
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
