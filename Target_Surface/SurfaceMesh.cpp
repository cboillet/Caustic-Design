#include "global.h"
#include <math.h>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*local*/
#include "SurfaceMesh.h"
#include "utils.h"



Mesh::Mesh(vector<Vertex> vertices, vector<GLuint> indices, vector<Texture> textures)
{
    Vertex vertex;
    glm::vec3 vect;
    this->vertices = vertices;
    this->indices = indices;
    this->textures = textures;
    nbMeshLayout = 2; //nb of triangle mesh in the space layout

    GLenum code = glewInit();
    if(code != GLEW_OK)
    {
        fprintf(stderr, "impossible d'initialiser GLEW : %s\n",
                        glewGetErrorString(code));
    }

    for (int i = 0; i<3; i++){
        vect.x = 0;
        vect.y = 0;
        vect.z = 0;
        vertex.Position = vect;

        vect.x = 0;
        vect.y = 0;
        vect.z = 0;
        vertex.Normal = vect;
        stopVertex.push_back(vertex);
    }

}

float Mesh::vectorNorm(Vertex v1, Vertex v2){
    int sum = pow((v2.Position.z-v1.Position.z),2)+pow((v2.Position.y-v1.Position.y),2)+pow((v2.Position.x-v1.Position.x),2);
    return sqrt(sum);
}

vector<int> Mesh::longerSegment(Vertex v0, Vertex v1, Vertex v2, int first){
    vector<int> ind;
    int max = vectorNorm(v0,v1);
    if(vectorNorm(v1,v2)>max) {
        ind.push_back(first + 1);
        ind.push_back(first + 2);
        ind.push_back(first); //no concerned by longest segment
    }
    else if(vectorNorm(v0,v2)>max) {
        ind.push_back(first);
        ind.push_back(first + 2);
        ind.push_back(first + 1);
    }
    else {
        ind.push_back(first);
        ind.push_back(first + 1);
        ind.push_back(first + 2);
    }
    return ind;
}

bool Mesh::adjacent(vector<Vertex> vec1, vector<Vertex> vec2){
    bool adj = false;
    int equal;
    for(int i=0; i<3; i++){
        for (int j=0; j<3; j++)
            if ((vec1[i].Position.x == vec2[j].Position.x) && (vec1[i].Position.y == vec2[j].Position.y) && (vec1[i].Position.z == vec2[j].Position.z)) equal++;
    }
    if (equal >= 2) adj = true;
    return adj;
}

bool Mesh::equal(vector<Vertex> vec1, vector<Vertex> vec2){
    bool adj = false;
    int equal;
    for(int i=0; i<3; i++){
        for (int j=0; j<3; j++)
            if ((vec1[i].Position.x == vec2[j].Position.x) && (vec1[i].Position.y == vec2[j].Position.y) && (vec1[i].Position.z == vec2[j].Position.z)) equal++;
    }
    if (equal >= 3) adj = true;
    return adj;
}

bool Mesh::compareArea(vector<Vertex> vec1, vector<Vertex> vec2){
    int a1,a2;
    vector<int> ia1 = longerSegment(vec1[0],vec1[1],vec1[2],0);
    vector<int> ia2 = longerSegment(vec2[0],vec2[1],vec2[2],0);
    int edgea1 = outTriplet(ia1,0,2);
    int edgea2 = outTriplet(ia2,0,2);
    a1 = vectorNorm(vec1[ia1[0]],vec1[edgea1])*vectorNorm(vec1[ia1[1]],vec1[edgea1])/2;
    a2 = vectorNorm(vec2[ia2[0]],vec2[edgea2])*vectorNorm(vec2[ia2[1]],vec2[edgea2])/2;
    if (a2>=a1) return true;
    else return false;
}


