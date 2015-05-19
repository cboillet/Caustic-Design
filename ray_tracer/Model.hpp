#ifndef MODEL_H
#define MODEL_H

#include "TriMesh.hpp"
#include "Image.hpp"

#include "glm/glm.hpp"
#define GLM_FORCE_RADIANS
#include "glm/gtx/unsigned_int.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/matrix_inverse.hpp"
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

class Model{

public:
    struct Material{
      // ambient color
      glm::vec4 ambient;
      // diffuse color
      glm::vec4 diffuse;
      // specular color
      glm::vec4 specular;
        // for ray tracer
      float transmit;
      // shininess
      float shininess;
    };


    Model();
    ~Model(){};

    glm::mat4* getModelMatrix();

    void loadMesh(std::string filename, TriMesh::PolygonWinding winding);
    void randomizeMaterial();

    void rotate(glm::mat4 newRotation);

    void translate(glm::vec3 newTranslation);

    void scale(float newScaleFactor);

    void createModelMatrix();

    void draw(void);

    void setSelected(bool isSelected){mesh.setSelected(isSelected);}
    bool getSelected(){return mesh.getSelected();}

    Material* getMaterial(){return &material;}

    Image* getImage(){return &texture;}

    TriMesh* getMesh(){return &mesh;}

    float getScaling(){return scaleFactor;}

private:

    void randomizeSingleMaterial(glm::vec4 & materialColor);

    TriMesh mesh;
    Material material;
    Image texture;

    glm::mat4 modelMatrix;
    glm::mat4 rotation;
    glm::vec3 translation;
    float scaleFactor;

    int modelIndex;

};

#endif //MODEL_H
