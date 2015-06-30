#ifndef SURFACEMODEL_H
#define SURFACEMODEL_H

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*openGL*/
#include <GL/glew.h>
#include <GL/gl.h>
/*local*/
#include "SurfaceMesh.h"

using namespace std;


GLint TextureFromFile(const char* path, string directory);

class Model{
    public:
        /*  Functions   */
        Model(GLchar* path);
        ~Model(){}
//        void Draw(Shader shader);
    private:
        /*  Model Data  */
        vector<Mesh> meshes;
        Mesh mesh;
        string directory;
        /*  Functions   */
        void loadModel(string path);
        void processNode(aiNode* node, const aiScene* scene);
        Mesh processMesh(aiMesh* mesh, const aiScene* scene);
        vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);
};

#endif // SURFACEMODEL_H
