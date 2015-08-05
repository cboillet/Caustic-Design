#ifndef SURFACEMESH_H
#define SURFACEMESH_H
#pragma comment(lib, \"assimp.lib\")

#include "global.h"
#include <math.h>
#include <vector>
#include <string>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
#include "utils.h"

using namespace std;

struct Vertex {
    glm::vec3 Position;
    glm::vec3 Normal;
    glm::vec2 TexCoords;
};

struct Texture {
    GLuint id;
    string type;
};


class Mesh {
public:
        /*  Mesh Data  */
        vector<Vertex> vertices;
        vector<Vertex> vertices2;
        vector<Vertex> verticestemp;
        vector<Vertex> stopVertex;
        vector<GLuint> indices;
        vector<Texture> textures;
        int nbMeshLayout;
       //int nbTriangle;
        /*  Render data  */
        GLuint VAO, VBO, EBO;
        /*  Functions  */

        Mesh(vector<Vertex> vertices, vector<GLuint> indices, vector<Texture> textures);
        Mesh(){}
        //void Draw(Shader shader);
        int getSpatialLayout(int nbFace);
        void generateTriangles();
        void setUpMesh(int nbvertices);
        float vectorNorm(Vertex v0, Vertex v1);
        vector<int> longerSegment(Vertex v0, Vertex v1, Vertex v2, int first);
        void recursiveInsert(vector<Vertex> vec, int nbTriangle, int up, int depth);
        void recursiveInsert2(vector<Vertex> vecteurList, int nbTriangle, int up, int depth);
        vector<Vertex> insertRecursiveVertex(vector<int> ind, Vertex v1,vector<Vertex> vec);
        void insertVertex(int& i1, int& i2, Vertex v1, vector<Vertex> vec);
        bool adjacent(vector<Vertex> vec1, vector<Vertex> vec2);
        bool equal(vector<Vertex> vec1, vector<Vertex> vec2);
        bool compareArea(vector<Vertex> vec1, vector<Vertex> vec2); //return true if vec2 bigger or equal than vec1
        void parseDiffMesh(vector<Vertex> v1, vector<Vertex> v2);
        //void findEqual(int[] v, vector<Vertex> v2, vector<Vertex>::iterator it);
};

class TargetSurfaceMesh : public Mesh {
    double m_dx; //corresponding to the one in the scene
    double m_dy;

public:
    TargetSurfaceMesh(double dx, double dy){}
    TargetSurfaceMesh(){}
    ~TargetSurfaceMesh(){}
    double get_dx() const { return m_dx; }
    double get_dy() const { return m_dy; }
};

#endif // SURFACEMESH_H
