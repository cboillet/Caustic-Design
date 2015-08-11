#ifndef UTILS_H
#define UTILS_H
#include"global.h"
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
#include <math.h>
#include <vector>
#include <string>

using namespace std;


int outTriplet(vector<int> vec, int begin, int end);
bool floatEquals(float val1, float val2);


//calculate the normals depending on x for the cost function

//void calculateFaceNormals(std::vector<glm::vec3> &normals, vector<glm::vec3> positions)
//{
//    for (uint i=0; i<indices.size(); i++)
//    {
//        glm::vec3 v1, v2;
//        v1 = positions[(i+1)%3] - positions[(i)%3];
//        v2 = positions[(i+2)%3] - positions[(i)%3];

//        normals.push_back(glm::normalize(glm::cross(v1, v2)));
//    }
//}

//vector<glm::vec3> calculateVertexNormals(vector<glm::vec3> positions)
//{
//    vector<glm::vec3> faceNormals;
//    vector<glm::vec3> normales

//    calculateFaceNormals(faceNormals, positions);

//    //aF contains indices of adjacend faces per vertex
//    vector<vector<unsigned int> > aF;
//    aF.resize(vertices.size());
//    for(uint i = 0; i < indices.size(); i++){
//        aF[indices[i][0]].push_back(i);
//        aF[indices[i][1]].push_back(i);
//        aF[indices[i][2]].push_back(i);
//    }

//    //interpolate normals of adjacend faces per vertex
//    for(uint i = 0; i< aF.size(); i++){

//        glm::vec3 vertexNormal = glm::vec3(0);
//        for(uint j= 0; j < aF[i].size(); j++){

//            // find out which vertex of the current face is the vertex we are currently looking at
//            // aF[i] is a list of faces (aka a list of indices of the indices-vector)
//            // so indices[aF[i][j]] is a glm::vec3 that contains one face
//            // and the current vertex is vertex[i]
//            int thisVertexIndex = -1;
//            for(int vIndex=0; vIndex < 3; vIndex++)
//            {
//                glm::vec3 pos1 = vertices[indices[aF[i][j]][vIndex]].Position;
//                glm::vec3 pos2 = vertices[i].Position;
//                if(glm::distance(pos1, pos2) < 0.0001f)
//                {
//                    thisVertexIndex = vIndex;
//                    break;
//                }

//            }

//            // we got index of our current vertex within the face, now get others
//            int other1 = (thisVertexIndex+1) % 3;
//            int other2 = (thisVertexIndex+2) % 3;

//            // create the vectors the represent the edges from current vertex to the other 2
//            glm::vec3 edge1 = vertices[indices[aF[i][j]][thisVertexIndex]].Position - vertices[indices[aF[i][j]][other1]].Position;
//            glm::vec3 edge2 = vertices[indices[aF[i][j]][thisVertexIndex]].Position - vertices[indices[aF[i][j]][other2]].Position;

//            // get angle between the edges
//            float incidentAngle = abs(glm::angle(glm::normalize(edge1), glm::normalize(edge2)));
//            if(incidentAngle > 180)
//               incidentAngle = 360 - incidentAngle;

//            // use that angle as weighting
//            vertexNormal += (faceNormals[aF[i][j]] * incidentAngle);
//        }
//        vertices[i].Normal = glm::normalize(vertexNormal);
//    }
//}


#endif // UTILS_H
