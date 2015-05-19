#include "Model.hpp"

using namespace std;

Model::Model(){
    mesh = TriMesh();
    rotation = glm::mat4(1);
    translation = glm::vec3(0);
    scaleFactor = 1.0f;
    randomizeMaterial();
}

void Model::loadMesh(string filename, TriMesh::PolygonWinding winding){
    mesh.setWinding(winding);
    mesh.loadOff(filename);
    mesh.center();
    mesh.unitize();
    mesh.calculateBoundingBox();
    mesh.calculateBoundingSphere();
    mesh.computeNormals();

    if(filename == "meshes/quad.off"){
        mesh.generateDefaultUV();
        texture.load("data/checker5.ppm");
    }else if(filename=="meshes/bunnysimple.off"){
        mesh.computeSphereUVs();
        texture.load("data/wood.ppm");
    }else if(filename == "meshes/killeroo.off"){
        mesh.computeSphereUVs();
        texture.load("data/flowers.ppm");
    }else{

        mesh.computeSphereUVs();
        if(!(texture.getName() == "data/world.ppm")){
            texture.load("data/world.ppm");
        }
    }
}

void Model::rotate(glm::mat4 newRotation){
    rotation = newRotation * rotation;
    createModelMatrix();
}

void Model::translate(glm::vec3 newTranslation){
    translation += newTranslation;
    createModelMatrix();
}

void Model::scale(float newScaleFactor){
    scaleFactor *= newScaleFactor;
    createModelMatrix();
}

void Model::randomizeMaterial(){
    srand(time(NULL));
    randomizeSingleMaterial(material.ambient);
    //randomizeSingleMaterial(material.diffuse);
    //randomizeSingleMaterial(material.specular);
    material.diffuse = glm::vec4(material.ambient);
    material.ambient = material.ambient * 0.2f;
    material.ambient[3] = 1.0f;
    material.diffuse = material.diffuse * 0.5f;
    material.diffuse[3] = 1.0f;

    material.specular = glm::vec4(1);
    material.specular = material.specular * 0.25f;
    material.specular[3] = 1.0f;

    material.transmit = 0.2f;
    material.shininess = rand() % 100;
}

void Model::randomizeSingleMaterial(glm::vec4 &materialColor){
    for(unsigned int i=3; i--;){
        materialColor[i] = (float)(rand()%255);
    }
    materialColor /= 255.f;
    materialColor[3] = 1.f;
}

void Model::createModelMatrix(){
    modelMatrix = glm::mat4(1);
    modelMatrix = glm::translate(modelMatrix, translation);
    modelMatrix *= rotation;
    modelMatrix = glm::scale(modelMatrix, glm::vec3(scaleFactor));
}

glm::mat4* Model::getModelMatrix() {
    return &modelMatrix;
}

void Model::draw(void){
    mesh.draw();
}

