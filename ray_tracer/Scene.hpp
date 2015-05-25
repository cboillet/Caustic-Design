#ifndef SCENE_HPP
#define SCENE_HPP

#include "TriMesh.hpp"
#include "Model.hpp"
#include "Image.hpp"

#include <glm/glm.hpp>
#include <string>


using namespace glm;
using namespace std;

namespace Scene{

    struct LightSource{
      //position in view space
      glm::vec4 position;
      // ambient color
      glm::vec4 ambient;
      // diffuse color
      glm::vec4 diffuse;
      // specular color
      glm::vec4 specular;
    };

    void render();
    glm::vec3 calculatePointColor(Ray & ray, glm::vec3 &intersectionPoint, glm::vec3 & baryIntersection, Model & model, int faceIndex, int recDepth);
    bool isInShadow(glm::vec3 & intersectionPoint);

    void loadModel(string path, TriMesh::PolygonWinding winding);

    void normalizeColor(glm::vec3 & color);


    extern vector<Model> viewElements;

    extern float lightStep;
    extern LightSource light;
    //extern Material material;

    extern int currentModelIndex;
    extern int currentElementIndex;

    // field of view (in degrees)
    extern GLfloat fov;
    // camera matrix
    extern mat4 cameraMatrix;
    //projection matrix
    extern mat4 projectionMatrix;

    // world rotation..
    extern glm::mat4 worldRotation;
    extern float worldXRotation;
    extern float worldYRotation;

    extern int recursionDepth;
    extern bool texturingEnabled;


    // for camera rotation
    extern float eyeXRot;// = 0.0f;
    extern float eyeYRot;// = 0.0f;
    extern glm::vec3 initialEyePos;// = vec3(0.0, 0.0, -8.0);
    extern glm::vec3 eye;// = vec3(initialEyePos);

    //render resolution
    extern unsigned int xRes;
    extern unsigned int yRes;

    //samples per pixel
    extern float samples;

    void setLighting();
    void setMaterial();

    void modifyLightColor(unsigned int colorIndex, float change);

    // TODO: this does not make sense for several objects. only for debug
    glm::mat4 getModelMatrix();
    glm::mat4 getModelMatrix(unsigned int index);

    void rotate(glm::mat4 newRotation);
    void translate(glm::vec3 newTranslation);
    void scale(float scaleFactor);

    void nextModel();
    void previousModel();
    TriMesh::PolygonWinding getModelWinding(std::string model);

    void updateModelMatrix();

    void addNewModel();
    void deleteCurrentElement();
    void nextElement();
    void previousElement();

    extern Image image;
}


#endif // SCENE_HPP
