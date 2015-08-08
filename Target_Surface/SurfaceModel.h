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
        Model();
        ~Model(){}
        void loadModel(string path, float scaling);
        void loadReceiverLightPoints(QString path);
//        void Draw(Shader shader);
        int findSurfaceMesh();
        void printAllVertices();
        void exportModel(std::string filename);
        vector<glm::vec3> getLightRayPositions(){ return receiverLightPositions; }
        float getFocalLength() { return focalLength; }
        void rescaleMeshes(float oldScale, float newScale);
        /*  Model Data  */
        vector<Mesh> meshes; //in our case just one Mesh -> this is for more complex models
        vector<glm::vec3> desiredNormals;
        vector<glm::vec3> currentNormals;

    protected:
        aiScene* scene;
        float focalLength;
        vector<glm::vec3> receiverLightPositions;

        Mesh SurfaceMesh;
        //Mesh mesh;
        string directory;
        /*  Functions   */
        virtual void loadToSurface(int index){}
        virtual void processNode(aiNode* node, const aiScene* scene, float scaling);
        virtual Mesh processMesh(aiMesh* mesh, const aiScene* scene, float scaling);
        virtual vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);

        void populateMesh(){}
        void compute(){}
        void clean(){
            meshes.clear();
        }
};

#endif // SURFACEMODEL_H
