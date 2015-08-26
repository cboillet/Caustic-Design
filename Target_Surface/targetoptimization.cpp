#include "targetoptimization.h"

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;





/**************************************************************/
/************** Old Versions of Cost-Functors *****************/
/**************************************************************/

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
        vertices = m->meshes[0].faceVertices;
        edgeVertices = m->meshes[0].faceVerticesEdge;
        edgeIndices = m->meshes[0].edgeIndices;
        edgeNoEdgeMapping = m->meshes[0].createEdgeToNoEdgeMapping();
        //noEdgeEdgeMapping = m->meshes[0].createNoEdgeToEdgeMapping();
        edgeAdjacentFaces = m->meshes[0].edgeAdjacentFaces;
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

        //std::cerr << "done " << std::endl;

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

            // get angle between the edges
            //float incidentAngle = abs(glm::angle(glm::normalize(edge1), glm::normalize(edge2)));
            //if(incidentAngle > 180)
            //   incidentAngle = 360 - incidentAngle;*/

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
    vector<Vertex*> vertices;
    vector<Vertex*> edgeVertices;
    vector<int> edgeNoEdgeMapping;
    vector<int> noEdgeEdgeMapping;
    vector<glm::uvec3> edgeIndices;
    vector<vector<uint> > edgeAdjacentFaces;
    float weight;

};


class CostFunctorEint4{
public:
    CostFunctorEint4(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {
    }

    bool operator()(const double* const vertex,
                                          const double* const neighbor1,
                                          const double* const neighbor2,
                                          const double* const neighbor3,
                                          const double* const neighbor4,
                                          const double* const neighbor5,
                                          const double* const neighbor6,
                                          double* residual) const
    {

        if(false)
        {
            residual[0] = 0;
            return true;
        }

        // -- preparation
        const double* neighbors[6];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;


        double** faceNormals = new double*[neighborMap.size()/2];


        // -- face normals
        for(uint i=0; i<neighborMap.size(); i+=2)
        {
            int faceIndex = i/2;
            faceNormals[faceIndex] = new double[3];
            double* faceNormal = faceNormals[faceIndex];

            const double* n1 = neighbors[neighborMap[i]];
            const double* n2 = neighbors[neighborMap[i+1]];
            calcFaceNormal(vertex, n1, n2, faceNormal);
        }

        // -- vertex normal
        double* vertexNormal = new double[3];
        calcVertexNormal(vertex, vertexNormal, faceNormals, neighbors);

        double vX = vertexNormal[0];
        double vY = vertexNormal[1];
        double vZ = vertexNormal[2];

        double dX = desiredNormal->x;
        double dY = desiredNormal->y;
        double dZ = desiredNormal->z;

        // evaluation
        double x = vertexNormal[0] - (desiredNormal->x);
        double y = vertexNormal[1] - (desiredNormal->y);
        double z = vertexNormal[2] - (desiredNormal->z);

        double res = sqrt(x*x + y*y + z*z);
        //std::cout << "res = " << res << std::endl;

        std::cout << "vertex-normal: " << vX << ", " << vY << ", " << vZ << std::endl;
        std::cout << "desired-normal: " << dX << ", " << dY << ", " << dZ << std::endl << std::endl;

        residual[0] = res;


        // -- clean
        delete[] vertexNormal;

        for(uint i=0; i<neighborMap.size(); i+=2)
        {
            delete[] faceNormals[i/2];
        }

        delete[] faceNormals;

        return true;
    }

    void calcFaceNormal(const double* const v1, const double* const v2, const double* const v3, double* result) const
    {
        double vertex1[3];
        double vertex2[3];

        vertex1[0] = v2[0] - v1[0];
        vertex1[1] = v2[1] - v1[1];
        vertex1[2] = v2[2] - v1[2];

        vertex2[0] = v3[0] - v1[0];
        vertex2[1] = v3[1] - v1[1];
        vertex2[2] = v3[2] - v1[2];

        cross(vertex1, vertex2, result);

        //std::cerr << "result = [" << result[0] << ", " << result[1] << ", " << result[2] << "]" << std::endl;
        //std::cerr << "foo" << std::endl;
    }

    void calcVertexNormal(const double* vertex, double* result, double** faceNormals, const double** neighbors) const
    {
        result[0] = result[1] = result[2] = 0;

        double edge1Res[3];
        double edge2Res[3];
        double edge1[3];
        double edge2[3];
        double edge1Sub[3];
        double edge2Sub[3];

        for(uint i=0; i<neighborMap.size(); i+=2)
        {

            for (uint j=0; j<3; j++)
            {
                edge1[j] = vertex[j];
                edge2[j] = vertex[j];
                edge1Sub[j] = neighbors[neighborMap[i]][j];
                edge2Sub[j] = neighbors[neighborMap[i+1]][j];
            }


            for(uint j=0; j<3; j++)
            {
                edge1Res[j] = edge1[j] - edge1Sub[j];
                edge2Res[j] = edge2[j] - edge2Sub[j];
            }

            normalize(edge1Res);
            normalize(edge2Res);

            double incidentAngle = angle(edge1Res, edge2Res);

            double* faceNormal = faceNormals[i/2];

            // use that angle as weighting
            for (uint j=0; j<3; j++)
            {
                double val = faceNormal[j];
                result[j] += val * incidentAngle;
            }

        }

        normalize(result);
    }


    double angle(double* v1, double* v2) const
    {
        return acos(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    }

    void normalize(double* v) const
    {
        double sum = 0;
        for (uint i=0; i<3; i++)
        {
            sum += v[i] * v[i];
        }

        sum = sqrt(sum);

        for (uint i=0; i<3; i++)
        {
            v[i] = v[i]/sum;
        }
    }


    void cross(double* v1, double* v2, double* result) const
    {
        result[0] = (v1[1]*v2[2] - v1[2]*v2[1]);
        result[1] = (v1[2]*v2[0] - v1[0]*v2[2]);
        result[2] = (v1[0]*v2[1] - v1[1]*v2[0]);
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

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



/********************************************************/
/************** TargetOptimization Class ****************/
/********************************************************/

TargetOptimization::TargetOptimization()
{
    //computeNormals.reserve(m.meshes[0].selectVerticesMeshFaceEdge().size());
}

TargetOptimization::~TargetOptimization(){
    //if (model != NULL) delete[] model;
}

void TargetOptimization::gatherVertexInformation(Vertex *vertex, uint vertexIndex, vector<int> &neighborList, vector<int> &neighborMap)
{
    Mesh * mesh = &model->meshes[0];
    vector<uint> adjacentFaces = mesh->edgeAdjacentFaces[vertexIndex];

    for(uint faceIndex = 0; faceIndex < adjacentFaces.size(); faceIndex++)
    {
        glm::uvec3 face = model->meshes[0].edgeIndices[adjacentFaces[faceIndex]];
        int thisVertexIndex = 0;
        for (int vIndex = 0; vIndex < 3; vIndex++)
        {
            Vertex * v = mesh->faceVerticesEdge[face[vIndex]];

            if(v == vertex){
                thisVertexIndex = vIndex;
                break;
            }
        }

        int other[] = {
            (thisVertexIndex+1)%3,
            (thisVertexIndex+2)%3
        };


        for(uint j=0; j<2; j++)
        {
            std::vector<int>::iterator it;
            it = std::find(neighborList.begin(), neighborList.end(), face[other[j]]);
            if(it != neighborList.end()) {
                // push back index of the neighbor
                neighborMap.push_back(it - neighborList.begin());
            } else {
                // push back neighbor and the index of it
                neighborList.push_back(face[other[j]]);
                neighborMap.push_back(neighborList.size()-1);
            }
        }
    }

}


void TargetOptimization::addResidualBlocks(Problem *problem, uint vertexIndex, vector<int> &neighbors, vector<int> &neighborMap, double *vertices)
{

    // EBar only calculates the distance from the vertex to the receiver-plane. So passing x-coordiate of receiving plane is sufficient
    glm::vec3 receiverPos = glm::vec3(model->meshes[0].getMaxX() + model->getFocalLength(), 0, 0);
    CostFunction* cost_function_ebar =
            new AutoDiffCostFunction<CostFunctorEbar2, 1, 3>(new CostFunctorEbar2(&receiverPos));

    problem->AddResidualBlock(
                cost_function_ebar,
                NULL, // no loss function
                &vertices[vertexIndex]
                );


    // EDir depends on the original position
    glm::vec3 pos; // TODO get original position
    CostFunction* cost_function_edir =
            new AutoDiffCostFunction<CostFunctorEdir2, 1, 3>(new CostFunctorEdir2(&pos));

    problem->AddResidualBlock(
                cost_function_edir,
                NULL,
                &vertices[vertexIndex]
                );


    // now we need to take care of the different amount of neighbors for EInt and EReg
}


void TargetOptimization::runTest(Renderer* renderer)
{

    Mesh * mesh = &model->meshes[0];


    vector<vector<int> > neighborsPerVertex;
    neighborsPerVertex.resize(mesh->faceVerticesEdge.size());

    vector<vector<int> > neighborMapPerVertex;
    neighborMapPerVertex.resize(mesh->faceVerticesEdge.size());

    bool allCaptured = true;
    for(uint i=0; i<mesh->faceVerticesEdge.size(); i++)
    {
        Vertex * v = mesh->faceVerticesEdge[i];
        vector<int> neighbors;
        vector<int> neighborMap;

        gatherVertexInformation(v, i, neighbors, neighborMap);

        uint nNeighbors = neighbors.size();
        if(nNeighbors != 2 && nNeighbors != 3 && nNeighbors != 4 && nNeighbors != 5 && nNeighbors != 6 && nNeighbors != 7 && nNeighbors != 8)
        {
            std::cout << "gathered information of vertex " << i << std::endl << "\tvertex = " << v->Position.x << ", " << v->Position.y << ", " << v->Position.z << std::endl;
            std::cout << "\tnNeighbors = " << neighbors.size() << std::endl;
            std::cout << "\tnFaces = " << (neighborMap.size()/2) << std::endl << std::endl;
            allCaptured = false;
        }

        neighborsPerVertex[i] = neighbors;
        neighborMapPerVertex[i] = neighborMap;
    }

    if(allCaptured)
        std::cout << "all possible neighbor-numbers captured" << std::endl;


    // prepare model and mesh
    mesh->calculateVertexNormals();
    model->computeLightDirectionsScreenSurface();
    model->fresnelMapping();

    // put all positions in one big list that we access later
    double* vertices = new double[3*mesh->faceVerticesEdge.size()];
    for(uint i=0; i<mesh->faceVerticesEdge.size(); i++)
    {
        glm::vec3 * pos = &mesh->faceVerticesEdge[i]->Position;

        vertices[3*i + 0] = pos->x;
        vertices[3*i + 1] = pos->y;
        vertices[3*i + 2] = pos->z;
    }

    Problem prob;

    // iterate over all vertices and add the corresponding residual blocks
    for(uint i=0; i<neighborsPerVertex.size(); i++)
    {
        addResidualBlocks(&prob, i, neighborsPerVertex[i], neighborMapPerVertex[i], vertices);
    }




    if(true)
        return;


    int vertexIndex = 14;
    Vertex * vertex = &model->meshes[0].vertices[vertexIndex];
    vector<uint> adjacentFaces = model->meshes[0].adjacentFaces[vertexIndex];
    vector<int> neighborMap;
    vector<int> neighborList;

    // we are now getting all neighbors. Each only once.
    // and at the same time we are creating a mapping. So we kind of create a adjacentfaces with vetrices.
    // so the neighborMap always has two entries that are (in combination with current vertex) one face.
    // and they are in correct winding already. So in optimizer we don't have to take care
    for(uint faceIndex = 0; faceIndex < adjacentFaces.size(); faceIndex++)
    {
        glm::uvec3 face = model->meshes[0].indices[adjacentFaces[faceIndex]];
        int thisVertexIndex = 0;
        for (int vIndex = 0; vIndex < 3; vIndex++)
        {
            Vertex * v = &model->meshes[0].vertices[face[vIndex]];

            if(v == vertex){
                thisVertexIndex = vIndex;
                break;
            }
        }

        int other[] = {
            (thisVertexIndex+1)%3,
            (thisVertexIndex+2)%3
        };


        for(uint j=0; j<2; j++)
        {
            std::vector<int>::iterator it;
            it = std::find(neighborList.begin(), neighborList.end(), face[other[j]]);
            if(it != neighborList.end()) {
                // push back index of the neighbor
                neighborMap.push_back(it - neighborList.begin());
            } else {
                // push back neighbor and the index of it
                neighborList.push_back(face[other[j]]);
                neighborMap.push_back(neighborList.size()-1);
            }
        }
    }


    /*    if(neighborList.size() != 6)
        {
            std::cout << "neighborlist for vertex " << i << " size = " << neighborList.size() << std::endl;
            std::cout << "vertex " << i << " = " << vertex->Position.x << ", " << vertex->Position.y << ", " << vertex->Position.z << std::endl;
        }
        if(neighborMap.size() != (2*neighborList.size()))
        {
            std::cout << "vertex " << i << "-- neighbors: " << neighborList.size() << " neighbormap.size: " << neighborMap.size() << std::endl;
        }
    }*/

    //renderer->setNeigbors(neighborList, neighborMap);

    int nNeighbors = neighborList.size();
    std::cout << "nNeighbors = " << nNeighbors << std::endl;

    //double* vertices = new double[3*(nNeighbors+1)];

    vertices[0] = vertex->Position.x;
    vertices[1] = vertex->Position.y;
    vertices[2] = vertex->Position.z;

    uint nVertices = 3*(nNeighbors+1);
    std::cout << "nVertices = " << nVertices << std::endl;
    for (uint i = 1; i <= nNeighbors; i++)
    {
        glm::vec3 pos = model->meshes[0].vertices[neighborList[i-1]].Position;

        uint at = 3*i;
        std::cout << "adding @" << at << ": " << pos.x << ", " << pos.y << ", " << pos.z << std::endl;
        vertices[at] = pos.x;
        vertices[at+1] = pos.y;
        vertices[at+2] = pos.z;

    }

    // prepare model
    model->meshes[0].calculateVertexNormals();
    model->computeLightDirectionsScreenSurface();
    model->fresnelMapping();
    //Mesh mesh = model->meshes[0];
    int mappedIndex = mesh->getIndexTargetSurface(&model->meshes[0].vertices[vertexIndex]);

    CostFunction* cost_fuction_eint =
            // one residual, 7 parameters
            new AutoDiffCostFunction<CostFunctorEint6Neighbors, 1, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint6Neighbors(&model->desiredNormals[mappedIndex], neighborMap));
    CostFunction* cost_function_ereg =
            new AutoDiffCostFunction<CostFunctorEreg6Neighbors, 1, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEreg6Neighbors(model, renderer, neighborList));

    glm::vec3 pos = glm::vec3(x_sources[vertexIndex]);
    CostFunction* cost_function_edir =
            new AutoDiffCostFunction<CostFunctorEdir2, 1, 3>(new CostFunctorEdir2(&pos));

    glm::vec3 receiverPos = glm::vec3(mesh->getMaxX() + model->getFocalLength(), 0, 0);
    CostFunction* cost_function_ebar =
            new AutoDiffCostFunction<CostFunctorEbar2, 1, 3>(new CostFunctorEbar2(&receiverPos));

    Problem problem;

    problem.AddResidualBlock(
            cost_fuction_eint,
            NULL,
            vertices, // vertex
            &vertices[3], // and the neighbors..
            &vertices[6],
            &vertices[9],
            &vertices[12],
            &vertices[15],
            &vertices[18]);

    problem.AddResidualBlock(
            cost_function_ereg,
            NULL,
            vertices, // vertex
            &vertices[3], // and the neighbors..
            &vertices[6],
            &vertices[9],
            &vertices[12],
            &vertices[15],
            &vertices[18]);


    problem.AddResidualBlock(
            cost_function_edir,
            NULL,
            vertices);

    problem.AddResidualBlock(
            cost_function_ebar,
            NULL,
            vertices);


    glm::vec3 *ps[neighborList.size()];
    for(uint i=0; i<neighborList.size(); i++)
    {
        ps[i] = new glm::vec3(x_sources[neighborList[i]]);
        int index = 3* (i+1);
        CostFunction* fkt =
                new AutoDiffCostFunction<CostFunctorEdir2, 1, 3>(new CostFunctorEdir2(ps[i]));

        problem.AddResidualBlock(
                    fkt,
                    NULL,
                    &vertices[index]);
    }


    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    //options.linear_solver_type = ceres::ITERATIVE_SCHUR; //large bundle adjustment problems
    options.max_num_iterations = 500;
    options.dense_linear_algebra_library_type = ceres::LAPACK;
    //options.visibility_clustering_type = ceres::SINGLE_LINKAGE;
    //options.preconditioner_type = ceres::CLUSTER_TRIDIAGONAL; // fast preconditioner
    string error;
    if(!options.IsValid(&error))
    {
        std::cout << "Options not valid: " << error << std::endl;
    }

    Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << std::endl;


    vertex->Position.x = vertices[0];
    vertex->Position.y = vertices[1];
    vertex->Position.z = vertices[2];

    for (uint i=1; i<=nNeighbors; i++)
    {

        Vertex* v = &model->meshes[0].vertices[neighborList[i-1]];
        v->Position.x = vertices[i*3 + 0];
        v->Position.y = vertices[i*3 + 1];
        v->Position.z = vertices[i*3 + 2];
    }

    for (uint i=0; i<neighborList.size(); i++)
    {
        delete ps[i];
    }

    delete[] vertices;

    model->meshes[0].calculateVertexNormals();

    std::cout << "vertex-position 14: " << vertex->Position.x << ", " << vertex->Position.y << ", " << vertex->Position.z << std::endl;
    std::cout << "vertex-normal 14: " << vertex->Normal.x << ", " << vertex->Normal.y << ", " << vertex->Normal.z << std::endl;

    renderer->repaint();
}


void TargetOptimization::runOptimization(Model* m, Renderer* renderer){

    if(true)
    {
        model = m;
        for(int i=0; i<m->meshes[0].vertices.size(); i++){
            glm::vec3 v = (m->meshes[0].vertices[i].Position);
            x_sources.push_back(v);
        }
        runTest(renderer);
        return;
    }


    int numberIteration = 0;
    model=m;
    for(int i=0; i<m->meshes[0].faceVertices.size(); i++){
        glm::vec3 v = (m->meshes[0].faceVertices[i]->Position);
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
    options.max_num_iterations = 50;
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

    for (int i=0; i<n; i++)
    {
        x[i]->Position.x = x1[i];
        x[i]->Position.y = x2[i];
        x[i]->Position.z = x3[i];
    }


    delete[] x1;
    delete[] initial_x1;
    delete[] x2;
    delete[] initial_x2;
    delete[] x3;
    delete[] initial_x3;


    model->meshes[0].calculateVertexNormals();
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


