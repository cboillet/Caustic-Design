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
        void loadModel(string path);
        void loadReceiverLightPoints(QString path);
//        void Draw(Shader shader);
        int findSurfaceMesh();
        void printAllVertices();
        void exportModel(std::string filename);
        vector<glm::vec3> getLightRayPositions(){ return receiverLightPositions; }
        float getFocalLength() { return focalLength; }
        void setFocalLength(float newLength);
        void rescaleMeshes(float newScale);
        void modifyMesh();
        void setNormals(bool edge);
        void shootTestRay(glm::highp_dvec3 & direction, glm::highp_dvec3 & redirect, glm::highp_dvec3 & endpoint);
        void computeLightDirectionsScreenSurface();
        void fresnelMapping(); //compute the position on the surface to for the desired normals
        /*  Model Data  */
        vector<Mesh> meshes; //in our case just one Mesh -> this is for more complex models
        vector<glm::vec3> desiredNormals;
        vector<glm::vec3> screenDirections;
        vector<glm::vec3> currentNormals;
        vector<glm::vec3> incidentNormals;
        vector<glm::vec3> receiverLightPositions;

        float surfaceSize;

    protected:
        aiScene* scene;
        float focalLength;


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
