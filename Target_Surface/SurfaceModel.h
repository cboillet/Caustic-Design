#ifndef SURFACEMODEL_H
#define SURFACEMODEL_H

#include "global.h"
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*local*/
#include "SurfaceMesh.h"



using namespace std;


GLint TextureFromFile(const char* path, string directory);

class Model {
    public:
        /*  Functions   */
        Model(GLchar* path);
        Model(){}
        ~Model(){}
        void loadModel(string path);
//        void Draw(Shader shader);
        int findSurfaceMesh();
        void printAllVertices();
        void exportModel(std::string filename);
        /*  Model Data  */
        vector<Mesh> meshes; //in our case just one Mesh -> this is for more complex models

    protected:
        aiScene* scene;

        Mesh SurfaceMesh;
        //Mesh mesh;
        string directory;
        /*  Functions   */
        virtual void loadToSurface(int index){}
        virtual void processNode(aiNode* node, const aiScene* scene);
        virtual Mesh processMesh(aiMesh* mesh, const aiScene* scene);
        virtual vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);

        void populateMesh(){}
        void compute(){}
        void clean(){
            meshes.clear();
        }
};

#endif // SURFACEMODEL_H
