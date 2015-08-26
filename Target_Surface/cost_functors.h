#ifndef COST_FUNCTORS_H
#define COST_FUNCTORS_H


#include "ceres/ceres.h"
#include "glog/logging.h"
#include "glm/glm.hpp"

#include "global.h"

/******************************************************/
/*************** Ceres-Helper-Methods *****************/
/******************************************************/

using namespace std;

template<typename T> void cross(T* v1, T* v2, T* result);
template<typename T> void calcFaceNormal(const T* const v1, const T* const v2, const T* const v3, T* result);
template<typename T> T angle(T* v1, T* v2);
template<typename T> void normalize(T* v);
template<typename T> T evaluateInt(const T* const vertex, const T** const neighbors, uint nNeighbors, const vector<int> & neighborMap);
template<typename T> void calcVertexNormal(const T* vertex, T* result, T** faceNormals, const T** neighbors, const vector<int> & neighborMap);
template<typename T> T evaluateReg(const T** const allVertices, const float* L, uint nVertices);

/******************************************************/
/*************** Ceres-Helper-Methods *****************/
/******************************************************/

// For EInt
template<typename T> void cross(T* v1, T* v2, T* result)
{
    result[0] = T(v1[1]*v2[2] - v1[2]*v2[1]);
    result[1] = T(v1[2]*v2[0] - v1[0]*v2[2]);
    result[2] = T(v1[0]*v2[1] - v1[1]*v2[0]);
}

template<typename T> void calcFaceNormal(const T* const v1, const T* const v2, const T* const v3, T* result)
{
    T vertex1[3];
    T vertex2[3];

    vertex1[0] = v2[0] - v1[0];
    vertex1[1] = v2[1] - v1[1];
    vertex1[2] = v2[2] - v1[2];

    vertex2[0] = v3[0] - v1[0];
    vertex2[1] = v3[1] - v1[1];
    vertex2[2] = v3[2] - v1[2];

    cross(vertex1, vertex2, result);
}

template<typename T> T angle(T* v1, T* v2)
{
    return ceres::acos(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

template<typename T> void normalize(T* v)
{
    T sum = T(0);
    for (uint i=0; i<3; i++)
    {
        sum += v[i] * v[i];
    }

    sum = ceres::sqrt(sum);

    for (uint i=0; i<3; i++)
    {
        v[i] = v[i]/sum;
    }
}

template<typename T> T evaluateInt(const T* const vertex, const T** neighbors, uint nNeighbors, const vector<int> & neighborMap, const glm::vec3 * desiredNormal)
{

    T** faceNormals = new T*[neighborMap.size()/2];

    // -- face normals
    for(uint i=0; i<neighborMap.size(); i+=2)
    {
        int faceIndex = i/2;
        faceNormals[faceIndex] = new T[3];
        T* faceNormal = faceNormals[faceIndex];
        const T* n1 = neighbors[neighborMap[i]];
        const T* n2 = neighbors[neighborMap[i+1]];

        calcFaceNormal(vertex, n1, n2, faceNormal);
    }

    // -- vertex normal
    T* vertexNormal = new T[3];
    calcVertexNormal(vertex, vertexNormal, faceNormals, neighbors, neighborMap);


    // evaluation
    T x = vertexNormal[0] - T(desiredNormal->x);
    T y = vertexNormal[1] - T(desiredNormal->y);
    T z = vertexNormal[2] - T(desiredNormal->z);

    T res = ceres::sqrt(x*x + y*y + z*z);

    // -- clean
    delete[] vertexNormal;

    for(uint i=0; i<neighborMap.size(); i+=2)
    {
        delete[] faceNormals[i/2];
    }

    delete[] faceNormals;

    // -- done

    return res;
}

template<typename T> void calcVertexNormal(const T* vertex, T* result, T** faceNormals, const T** neighbors, const vector<int> & neighborMap)
{
    result[0] = result[1] = result[2] = T(0);

    T edge1Res[3];
    T edge2Res[3];
    T edge1[3];
    T edge2[3];
    T edge1Sub[3];
    T edge2Sub[3];

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

        T incidentAngle = angle(edge1Res, edge2Res);

        T* faceNormal = faceNormals[i/2];

        // use that angle as weighting
        for (uint j=0; j<3; j++)
        {
            T val = faceNormal[j];
            result[j] += val * incidentAngle;
        }

    }

    normalize(result);
}


// For EReg
template<typename T> T evaluateReg(const T** const allVertices, const float* L, uint nVertices)
{
    T res[3] = {T(0), T(0), T(0)};

    for(int i=0; i<nVertices; i++){
        for (int j=0; j<3; j++){
            res[j] += T(L[i]) * allVertices[i][j];
        }
    }

    return T(EREG_WEIGHT) * (ceres::abs(res[0]) + ceres::abs(res[1]) + ceres::abs(res[2]));

}


/**********************************************/
/************** Cost-Functors *****************/
/**********************************************/

/********* EInt *********/

class CostFunctorEint8Neighbors{
public:
    CostFunctorEint8Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          const T* const neighbor7,
                                          const T* const neighbor8,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[8];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;
        neighbors[6] = neighbor7;
        neighbors[7] = neighbor8;


        residual[0] = evaluateInt(vertex, neighbors, 8, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint7Neighbors{
public:
    CostFunctorEint7Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          const T* const neighbor7,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[7];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;
        neighbors[6] = neighbor7;


        residual[0] = evaluateInt(vertex, neighbors, 7, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint6Neighbors{
public:
    CostFunctorEint6Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[6];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;


        residual[0] = evaluateInt(vertex, neighbors, 6, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};

class CostFunctorEint5Neighbors{
public:
    CostFunctorEint5Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[5];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;


        residual[0] = evaluateInt(vertex, neighbors, 5, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint4Neighbors{
public:
    CostFunctorEint4Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[4];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;


        residual[0] = evaluateInt(vertex, neighbors, 4, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint3Neighbors{
public:
    CostFunctorEint3Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[3];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;


        residual[0] = evaluateInt(vertex, neighbors, 3, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint2Neighbors{
public:
    CostFunctorEint2Neighbors(glm::vec3 * desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[3];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;


        residual[0] = evaluateInt(vertex, neighbors, 3, neighborMap, desiredNormal);

        return true;
    }

private:
    glm::vec3 * desiredNormal;
    vector<int> neighborMap;

};

/********* EReg *********/

class CostFunctorEreg8Neighbors{
public:
    CostFunctorEreg8Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg8Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          const T* const neighbor7,
                                          const T* const neighbor8,
                                          T* residual) const
    {

        const T* allVertices[9];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;
        allVertices[4] = neighbor4;
        allVertices[5] = neighbor5;
        allVertices[6] = neighbor6;
        allVertices[7] = neighbor7;
        allVertices[8] = neighbor8;


        residual[0] = evaluateReg(allVertices, L, 9);

        return true;
    }

private:
    Model* model;
    float* L;
};

class CostFunctorEreg7Neighbors{
public:
    CostFunctorEreg7Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg7Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          const T* const neighbor7,
                                          T* residual) const
    {

        const T* allVertices[8];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;
        allVertices[4] = neighbor4;
        allVertices[5] = neighbor5;
        allVertices[6] = neighbor6;
        allVertices[7] = neighbor7;


        residual[0] = evaluateReg(allVertices, L, 8);

        return true;
    }

private:
    Model* model;
    float* L;
};


class CostFunctorEreg6Neighbors{
public:
    CostFunctorEreg6Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg6Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          T* residual) const
    {

        const T* allVertices[7];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;
        allVertices[4] = neighbor4;
        allVertices[5] = neighbor5;
        allVertices[6] = neighbor6;


        residual[0] = evaluateReg(allVertices, L, 7);

        return true;
    }

private:
    Model* model;
    float* L;
};

class CostFunctorEreg5Neighbors{
public:
    CostFunctorEreg5Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg5Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          T* residual) const
    {

        const T* allVertices[6];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;
        allVertices[4] = neighbor4;
        allVertices[5] = neighbor5;


        residual[0] = evaluateReg(allVertices, L, 6);

        return true;
    }

private:
    Model* model;
    float* L;
};


class CostFunctorEreg4Neighbors{
public:
    CostFunctorEreg4Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg4Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          T* residual) const
    {

        const T* allVertices[5];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;
        allVertices[4] = neighbor4;


        residual[0] = evaluateReg(allVertices, L, 5);

        return true;
    }

private:
    Model* model;
    float* L;
};


class CostFunctorEreg3Neighbors{
public:
    CostFunctorEreg3Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg3Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          T* residual) const
    {

        const T* allVertices[4];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;


        residual[0] = evaluateReg(allVertices, L, 4);

        return true;
    }

private:
    Model* model;
    float* L;
};


class CostFunctorEreg2Neighbors{
public:
    CostFunctorEreg2Neighbors(Model* m, Renderer* renderer, vector<int> neighbors): model(m){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }
    }

    ~CostFunctorEreg2Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          T* residual) const
    {

        const T* allVertices[3];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;


        residual[0] = evaluateReg(allVertices, L, 3);

        return true;
    }

private:
    Model* model;
    float* L;
};

/********* EBar *********/

class CostFunctorEbar2 {
public:
    CostFunctorEbar2 (glm::vec3* receiverPosition): receiverPosition(receiverPosition)
    {}

    template <typename T> bool operator()(const T* const vertex, T* e) const{

        T dth = T(EBAR_DETH);

        T distance = vertex[0] - T(receiverPosition->x);
        e[0] = T(EBAR_WEIGHT) * ceres::max(T(0), - ceres::log( (T(1)-distance) + dth) );

        return true;
    }

private:
    glm::vec3* receiverPosition;

};



/********* EDir *********/

class CostFunctorEdir2{
public:
    CostFunctorEdir2(glm::vec3 * s): source(s){
    }

    template <typename T> bool operator()(const T* const vertex, T* e) const{
        // x - proj(..) will only return the difference in y- and z-coordinates. we may ignore the rest
        T y = vertex[1] - T(source->y);
        T z = vertex[2] - T(source->z);
        e[0] = T(EDIR_WEIGHT) * (ceres::abs(y)+ceres::abs(z));
        return true;
    }

private:
    glm::vec3 * source;
    float weight;
};


#endif // COST_FUNCTORS_H
