#include "Scene.hpp"
#include "TriMesh.hpp"
#include "Assets.hpp"
#include "Context.hpp"
#include "Model.hpp"
#include "Ray.hpp"
#include "View3D.hpp"
#include "ViewRender.hpp"

#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_access.hpp>

#include <omp.h>

// the model
//TriMesh Scene::model;
vector<Model> Scene::viewElements = vector<Model>();

// field of view for camera
GLfloat Scene::fov = 60.0;
// camera position
mat4 Scene::cameraMatrix = glm::lookAt(vec3(0.0, 0.0, -8.0), vec3(0), vec3(0.0, 1.0, 0.0));
mat4 Scene::projectionMatrix;

// index of model being displayed
int Scene::currentModelIndex = 0;
int Scene::currentElementIndex = 0;

Scene::LightSource Scene::light;
//Scene::Material Scene::material;

static float lightDistance = 3.0;
float Scene::lightStep = 0.05;

glm::mat4 Scene::worldRotation = mat4(1);
float Scene::worldXRotation = 0.f;
float Scene::worldYRotation = 0.f;

float Scene::eyeXRot = 0.0f;
float Scene::eyeYRot = 0.0f;
glm::vec3 Scene::initialEyePos = vec3(0.0, 0.0, -8.0);
glm::vec3 Scene::eye = vec3(initialEyePos);

//TODO do we want to use window size here?
unsigned int Scene::xRes = 100;
unsigned int Scene::yRes = 100;

bool Scene::texturingEnabled = true;
int Scene::recursionDepth = 2;

//samples per pixel
float Scene::samples = 1.0f;

Image Scene::image(xRes, yRes);

void Scene::setLighting(){
    light.position= vec4(lightDistance, lightDistance, -lightDistance,1);
    light.ambient= vec4(1,1,1,1);
    light.diffuse= vec4(1,1,1,1);
    light.specular= vec4(0,0,0,0);
}

/*void Scene::setMaterial(){
    material.ambient= vec4(0.2,0.3,0.8,1);
    material.diffuse= vec4(0.2,0.3,0.8,1);
    material.specular= vec4(1, 1, 1, 1);
    material.shininess= 50;
}*/

void Scene::addNewModel(){
    if(viewElements.size() > 0){
        viewElements[currentElementIndex].setSelected(false);
    }

    // TODO: check what model-winding is!
    Model newModel = Model();
    newModel.setSelected(true);
    viewElements.push_back(newModel);
    currentElementIndex = viewElements.size()-1;
    loadModel(Assets::modelPath + Assets::models[currentModelIndex], getModelWinding(Assets::models[currentModelIndex]));
}

void Scene::loadModel(string path, TriMesh::PolygonWinding winding){
    if(viewElements.empty()) return;

    cout << "loading "<< path <<endl;
    viewElements[currentElementIndex].loadMesh(path, winding);
}

void Scene::nextModel() {
    currentModelIndex = (currentModelIndex+1)%Assets::modelAmount;
    loadModel(Assets::modelPath + Assets::models[currentModelIndex], getModelWinding(Assets::models[currentModelIndex]));
}

void Scene::previousModel(){
    if(currentModelIndex == 0){
        currentModelIndex = Assets::modelAmount-1;
    }else{
        currentModelIndex--;
    }
    loadModel(Assets::modelPath + Assets::models[currentModelIndex], getModelWinding(Assets::models[currentModelIndex]));
}

TriMesh::PolygonWinding Scene::getModelWinding(string model){
    if(model == "bunny.off" || model == "cow.off" || model == "cone.off" || model == "plane4x4.off" || model == "sphere.off"){
        return TriMesh::CCW;
    }else{
        return TriMesh::CW;
    }
}

void Scene::deleteCurrentElement(){
    if(viewElements.size() == 1){
        cerr << "must have at least one element left" << endl;
        return;
    }
    // remove element
    viewElements.erase(viewElements.begin()+currentElementIndex);

    currentElementIndex = 0;
    viewElements[currentElementIndex].setSelected(true);
}

void Scene::nextElement(){
    if(viewElements.empty()) return;

    viewElements[currentElementIndex].setSelected(false);

    currentElementIndex = (currentElementIndex+1)%viewElements.size();

    viewElements[currentElementIndex].setSelected(true);
}

void Scene::previousElement(){
    if(viewElements.empty()) return;

    viewElements[currentElementIndex].setSelected(false);

    if(currentElementIndex == 0){
        currentElementIndex = viewElements.size()-1;
    }else{
        currentElementIndex--;
    }

    viewElements[currentElementIndex].setSelected(true);
}

void Scene::updateModelMatrix(){
    //modelMatrix = translate(mat4(1), translation);
    //modelMatrix *= rotation;
}

mat4 Scene::getModelMatrix(){
    return (*viewElements[currentElementIndex].getModelMatrix());
}

mat4 Scene::getModelMatrix(unsigned int index){
    return (*viewElements[index].getModelMatrix());
}

void Scene::rotate(mat4 newRotation){
    viewElements[currentElementIndex].rotate(newRotation);
}

void Scene::translate(vec3 newTranslation){
    viewElements[currentElementIndex].translate(newTranslation);
}

void Scene::scale(float scaleFactor){
    viewElements[currentElementIndex].scale(scaleFactor);
}

void Scene::render(){

    // clean up previous (if any)
    ViewRender::removeOldPoints();

    // set size to current window dimensions
    xRes = ViewRender::screen.x*1.5f;
    yRes = ViewRender::screen.y*1.5f;
    image.setSize(xRes, yRes);

    cout << "\nstarted rendering..."<< endl;
    // get viewport
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);

    //generate primary rays
    unsigned int totalSamples = samples*xRes*yRes;
    cout << "reserving space for " << totalSamples << " rays, with " << sizeof(Ray) << " bytes/Ray " << endl;
    Ray* primaryRays = new Ray[totalSamples];

    float step = 1.0/samples;

    unsigned int pixel = 0;
    unsigned int x = 0;
    unsigned int y = 0;

    for(unsigned int i = totalSamples; i--;){
        //generate Ray
        //cout << "generate ray for pixel "<< int(i*step)<< endl;
        pixel = int(i*step);
        x = pixel%xRes;
        y = pixel/xRes;

        glm::vec3 nearPlanePos = glm::unProject(glm::vec3(((float)x/(float)xRes) * View3D::screen.x, ((float)y/(float)yRes) * View3D::screen.y, 0), cameraMatrix, projectionMatrix, glm::vec4(viewport[0],viewport[1],viewport[2], viewport[3]));
        glm::vec3 farPlanePos = glm::unProject(glm::vec3(((float)x/(float)xRes) * View3D::screen.x, ((float)y/(float)yRes) * View3D::screen.y, 1), cameraMatrix, projectionMatrix, glm::vec4(viewport[0],viewport[1],viewport[2], viewport[3]));

        primaryRays[int(i*step)] = Ray(nearPlanePos, normalize(farPlanePos-nearPlanePos));

    }

    cout << "rays generated.. calculating intersection" << endl;
    //trace the rays here

    vec3 intersection;
    vec3 baryPoint;


    vec3 cameraPosition = vec3(column(cameraMatrix, 3));
    int percentage = 0;
    int lastPercentage = -1;
    unsigned int done = 0;
    #ifndef WITH_OPENMP
    for(unsigned int i = totalSamples; i--;){
    #else
    #pragma omp parallel for shared(image, lastPercentage, percentage, done)
    for(unsigned int i = 0; i < totalSamples; i++){
    #endif
        percentage = int(100.0d/double(totalSamples)*(done));
        if(percentage > lastPercentage){
            lastPercentage = percentage;
            cout << percentage << " % "<< endl;
        }
        bool foundForRay = false;
        float bestDistance;
        vec3 bestFound;
        vec3 bestBary;
        int bestIndex;
        int meshIndex;
        int bestMeshIndex;
        for(unsigned j = viewElements.size(); j--;){
            mat4 modelMatrix = Scene::worldRotation * mat4(*viewElements[j].getModelMatrix());

            if((viewElements[j].getMesh())->intersectTriangle(primaryRays[i], modelMatrix, cameraPosition, intersection, baryPoint, meshIndex)){

                float distance = glm::distance(cameraPosition, intersection);

                if(foundForRay){
                    if(distance < bestDistance){
                        bestIndex = j;
                        bestDistance = distance;
                        bestFound = vec3(intersection);
                        bestBary = vec3(baryPoint);
                        bestMeshIndex = meshIndex;
                    }
                }else{
                    bestIndex = j;
                    bestDistance = distance;
                    bestFound = vec3(intersection);
                    bestBary = vec3(baryPoint);
                    bestMeshIndex = meshIndex;
                    foundForRay = true;
                }
            }
        }

        if(foundForRay){
            vec3 color = calculatePointColor(primaryRays[i], bestFound, bestBary, viewElements[bestIndex], bestMeshIndex, Scene::recursionDepth);
            ViewRender::addPoint(bestFound, color);
            image.setPixel(int(i/samples), color);
        }
        done++;
    }

    image.normalize();
    ViewRender::normalize();

    ViewRender::updateImage();
    Context::displayRenderResult();

    #ifdef WITH_OPENIMAGEIO
    image.saveToDisk("test.jpg");
    #else
    image.savePPM();
    #endif

    // after execution, free memory
    delete[] primaryRays;
    cout << "cleaned up. done rendering" << endl;

}

vec3 Scene::calculatePointColor(Ray &ray, vec3 &intersectionPoint, vec3 & baryIntersection, Model &model, int faceIndex, int recDepth){

    //vec3 cameraPosition = vec3(column(cameraMatrix, 3));
    mat4 modelMatrix = Scene::worldRotation * (*model.getModelMatrix());
    mat3 normalMatrix = mat3(transpose(inverse(cameraMatrix*modelMatrix)));
    mat3 normalModelMatrix = mat3(transpose(inverse(modelMatrix)));
    // for details refer to blinnPhong, note that if index == -1, the baryintersection holds the normal in world-coordinates

    // triNormal holds the normal in model-coordinates
    vec3 triNormal = faceIndex == -1 ? vec3(glm::inverse(modelMatrix) * vec4(baryIntersection,1)) : model.getMesh()->computeTriangleNormal(baryIntersection, faceIndex);
    // for diffuse light, transform model as in blinnphong
    vec3 triangleNormal = normalMatrix* triNormal;
    vec3 lightToPoint = vec3(cameraMatrix*Scene::light.position - cameraMatrix*vec4(intersectionPoint,1));

    float lightCos = glm::clamp(glm::dot(glm::normalize(triangleNormal), glm::normalize(lightToPoint)), 0.f, 1.f);
    vec3 diffuse;// = vec3((model.getMaterial()->diffuse * Scene::light.diffuse)) * lightCos;
    vec3 ambient;// = vec3(model.getMaterial()->ambient * Scene::light.ambient);

    if(!Scene::texturingEnabled){
        diffuse = vec3((model.getMaterial()->diffuse * Scene::light.diffuse)) * lightCos;
        ambient = vec3(model.getMaterial()->ambient * Scene::light.ambient);
    }else{
        vec2 uv;
        model.getMesh()->getUVCoords(baryIntersection, faceIndex, uv);

        vec4 color = model.getImage()->getRelativePixel(uv);

        diffuse = vec3(color * Scene::light.diffuse) * lightCos;
        ambient = vec3(color * Scene::light.ambient);
    }

    if(recDepth > 0){
        vec3 rayDir;
        ray.getDirection(rayDir);
        // here we don't want to be in camera space, but world space. so don't include camera in normal matrix
        Ray newRay = Ray(vec3(intersectionPoint), glm::normalize(glm::reflect(glm::normalize(rayDir),glm::normalize(normalModelMatrix * triNormal))));
        bool foundForRay = false;
        float bestDistance;
        vec3 bestFound;
        vec3 bestBary;
        int bestIndex;
        int meshIndex;
        int bestMeshIndex;
        vec3 intersection;
        vec3 baryPoint;
        for(unsigned j = viewElements.size(); j--;){
            mat4 newModelMatrix = Scene::worldRotation * mat4(*viewElements[j].getModelMatrix());
            // ray may be modified, so copy it
            Ray rayCopy = Ray(newRay);
            if((viewElements[j].getMesh())->intersectTriangle(rayCopy, newModelMatrix, intersectionPoint, intersection, baryPoint, meshIndex)){

                float distance = glm::distance(intersectionPoint, intersection);

                if(foundForRay){
                    if(distance < bestDistance){
                        bestIndex = int(j);
                        bestDistance = float(distance);
                        bestFound = vec3(intersection);
                        bestBary = vec3(baryPoint);
                        bestMeshIndex = int(meshIndex);
                    }
                }else{
                    bestIndex = int(j);
                    bestDistance = float(distance);
                    bestFound = vec3(intersection);
                    bestBary = vec3(baryPoint);
                    bestMeshIndex = int(meshIndex);
                    foundForRay = true;
                }
            }
        }

        if(foundForRay){

            vec3 color;
            // we have ambient + specular light in any case
            color = ambient + vec3(model.getMaterial()->specular) *calculatePointColor(newRay, bestFound, bestBary, viewElements[bestIndex], bestMeshIndex, recDepth-1);
            // we only have diffuse light if not in shadow
            if(!isInShadow(intersectionPoint)){
                color = vec3(color + diffuse);
            }


            return color;

        }else{
            vec3 color = ambient;
            if(!isInShadow(intersectionPoint)){
                color = vec3(color + diffuse);
            }



            //normalizeColor(color);
            return color;
        }

    }else{
        vec3 color = ambient;
        if(!isInShadow(intersectionPoint)){
             color = vec3(color + diffuse);
        }

        //normalizeColor(color);
        return color;
    }
}

bool Scene::isInShadow(vec3 & intersectionPoint){

    //return false;
    // we want to have intersectionPoint in world coordinates!

    vec3 pointToLight = glm::normalize(vec3(Scene::light.position) - intersectionPoint) ;
    Ray ray = Ray(intersectionPoint, pointToLight);

    // unused but needed for method
    vec3 intersection;
    vec3 baryPoint;
    int meshIndex;

    for(int j=0; j<viewElements.size(); j++){
        mat4 modelMatrix = Scene::worldRotation * (*viewElements[j].getModelMatrix());
        // ray may be modified in method
        Ray rayCopy = Ray(ray);
        if((viewElements[j].getMesh())->intersectTriangle(rayCopy, modelMatrix, intersectionPoint, intersection, baryPoint, meshIndex)){
            return true;
        }
    }

    return false;
}

void Scene::normalizeColor(vec3 & color){
    float max = 0.f;
    for(int i=0; i<3; i++)
        if(color[i] > max)
            max = color[i];
    if(max > 1.0f){
        for(int i=0; i<3; i++){
            color[i] /= max;
        }
    }
}

void Scene::modifyLightColor(unsigned int colorIndex, float step){

    vec4 min(0);
    vec4 max(1);

    light.ambient[colorIndex] += step;
    light.diffuse[colorIndex] += step;
    light.specular[colorIndex] += step;

    light.ambient = glm::clamp(light.ambient, min, max);
    light.diffuse = glm::clamp(light.diffuse, min, max);
    light.specular = glm::clamp(light.specular, min, max);
}
