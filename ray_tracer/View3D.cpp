#include "View3D.hpp"
#include "Scene.hpp"
#include "Assets.hpp"
#include "TriMesh.hpp"
#include "GLSLShader.hpp"
#include "Context.hpp"
#include "Model.hpp"

//convinience imports
#define GLM_FORCE_RADIANS

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <iostream>

#ifdef __APPLE__
  #include <OpenGL/gl.h>
#include <GL/freeglut.h>
#elif _WIN32
  #include "win32/glew.h"
#include "win32/freeglut.h"
#else
  #include <GL/glew.h>
#include <GL/freeglut.h>
#endif

using namespace glm;
using namespace std;


#define PI  3.14159265358979323846264338327950288f
#define RADIANS(x) (((x)*PI)/180.0f)

// current state of mouse action
static enum{
    ROTATE, SHIFT_XY, SHIFT_Z, SCALE, NO_DRAG, DRAW, ERASE, SCALEBRUSH, SMOOTHBRUSH, ROTATE_CAMERA
} drag= NO_DRAG;

static vec2 subwindow;

static GLSLShader shader;

static GLSLShader passpartoutShader;

static TriMesh passpartout("meshes/quad.off", false);

// the previous mouse position -- important for drag
vec2 previousMouse = vec2();
// screen dimensions
vec2 View3D::screen;

void View3D::loadShader(){
    cout << "in load shaders"<< endl;
    shader.loadVertexShader("shaders/preview.vert");
    shader.loadFragmentShader("shaders/blinnPhongReflection");
    shader.loadFragmentShader("shaders/preview.frag");
    shader.bindVertexAttrib("position", TriMesh::attribVertex);
    shader.bindVertexAttrib("normal", TriMesh::attribNormal);
    shader.link();

    passpartoutShader.loadVertexShader("shaders/passpartout.vert");
    passpartoutShader.loadFragmentShader("shaders/passpartout.frag");
    passpartoutShader.bindVertexAttrib("position", TriMesh::attribVertex);
    passpartoutShader.link();
}

void View3D::display(){

    glClearColor(0.5, 0.5, 0.5, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    shader.bind();

    // TODO: maybe one shader (and image?) per model?
    for(unsigned int i=0; i<Scene::viewElements.size(); i++){

        Model::Material* material = Scene::viewElements[i].getMaterial();
        mat4 modelMatrix = Scene::worldRotation * Scene::getModelMatrix(i);
        shader.setUniform("modelViewProjectionMatrix", Scene::projectionMatrix*Scene::cameraMatrix*modelMatrix);
        shader.setUniform("modelViewMatrix", Scene::cameraMatrix*modelMatrix);
        shader.setUniform("normalMatrix", mat3(transpose(inverse(Scene::cameraMatrix*modelMatrix))));

        shader.setUniform("lightSource.position", Scene::cameraMatrix * Scene::light.position);
        shader.setUniform("lightSource.ambient", Scene::light.ambient);
        shader.setUniform("lightSource.diffuse", Scene::light.diffuse);
        shader.setUniform("lightSource.specular", Scene::light.specular);

        shader.setUniform("material.ambient", material->ambient);
        shader.setUniform("material.diffuse", material->diffuse);
        shader.setUniform("material.specular", material->specular);
        shader.setUniform("material.shininess", material->shininess);

        Scene::viewElements[i].draw();
    }


    shader.unbind();

    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);

    passpartoutShader.bind();

    passpartoutShader.setUniform("renderSize", vec2(Scene::xRes, Scene::yRes));
    passpartoutShader.setUniform("screenSize", screen);
    passpartoutShader.setUniform("opacity", 0.6f);

    passpartout.draw();

    passpartoutShader.unbind();

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);



    glutSwapBuffers();

}

void View3D::reshape(int width, int height){
    subwindow = vec2(width, height);
    glViewport(0.0, 0.0, subwindow.x, subwindow.y);


    float ratio = (float)width/(float)height;
    Scene::projectionMatrix = perspective(radians(Scene::fov), ratio, 1.0f, 256.0f);

    screen = vec2(width, height);

}

void View3D::keyPressed(unsigned char key, int x, int y){
    (void)x;
    (void)y;

    switch (key){
        case 'q':
        case 'Q':
            exit(EXIT_SUCCESS);

        case 'm':
            Scene::nextModel();
            break;
        case 'M':
            Scene::previousModel();
            break;

        case 's':
            // TODO: this is only for debugging
            //Scene::viewElements.setSelected(!Scene::viewElements.getSelected());
            break;

        case 'e':
            Scene::nextElement();
            break;
        case 'E':
            Scene::previousElement();
            break;
        case 'd':
            Scene::deleteCurrentElement();
            break;
        case 'n':
            Scene::addNewModel();
            break;
        case 'r':
        case 'R':
            resetCamera();
        break;
        case 't':
            Scene::modifyLightColor(0, Scene::lightStep);
            break;
        case 'T':
            Scene::modifyLightColor(0, -Scene::lightStep);
            break;
        case 'g':
            Scene::modifyLightColor(1, Scene::lightStep);
            break;
        case 'G':
            Scene::modifyLightColor(1, -Scene::lightStep);
            break;
        case 'b':
            Scene::modifyLightColor(2, Scene::lightStep);
            break;
        case 'B':
            Scene::modifyLightColor(2, -Scene::lightStep);
            break;
    }

    Context::display3DView();
}


void View3D::mouseDragged(int x, int y){

    glm::mat4 rotX;
    glm::mat4 rotY;

    vec2 v= (vec2(x,y) - previousMouse) / screen;

    switch(drag){
        case ROTATE:
            if(length(v)==0) break;
            Scene::rotate(glm::rotate(mat4(1), RADIANS(180 * length(v)), normalize(vec3(-v.y, v.x, 0))));
            break;
        case SHIFT_XY:
            //Scene::translation.x -= 3.3*v.x;
            //Scene::translation.y -= 3.3*v.y;
            Scene::translate (-vec3(3.3*v.x, 3.3*v.y, 0));
            break;
        case SHIFT_Z:
            Scene::translate(vec3(0.0, 0.0, 3.3*sign(dot(v, vec2(1,1))) * length(v)));
            break;
        case SCALE:
            Scene::scale(1.0+10*v.y);
            break;
        case ROTATE_CAMERA:
            /*Scene::eyeYRot -= v.x*3.f;
            Scene::eyeXRot += v.y*3.f;
            rotX = glm::rotate(mat4(1), Scene::eyeYRot, vec3(0,1,0));
            rotY = glm::rotate(mat4(1), Scene::eyeXRot, vec3(1,0,0));
            Scene::eye = vec3(rotX * rotY * vec4(Scene::initialEyePos,1));
            updateCamera();*/
            Scene::worldXRotation -= v.y*3.f;
            Scene::worldYRotation += v.x*3.f;
            rotY = glm::rotate(mat4(1), Scene::worldYRotation, vec3(0,1,0));
            rotX = glm::rotate(mat4(1), Scene::worldXRotation, vec3(1,0,0));
            Scene::worldRotation = rotX *rotY;
            break;
        default:
            break;
    }

    previousMouse= vec2(x, y);

    Scene::updateModelMatrix();
    Context::display3DView();
}

void View3D::mousePressed(int button, int state, int x, int y){

    int modifier;

    switch(button){
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN){
        previousMouse= vec2(x, y);
        modifier = glutGetModifiers();
        if(modifier & GLUT_ACTIVE_CTRL)
            drag = SHIFT_XY;
        else if(modifier & GLUT_ACTIVE_SHIFT)
            drag = SHIFT_Z;
        else
            drag = ROTATE;
        }
        else if(state == GLUT_UP){
            drag = NO_DRAG;
        }
    break;
    case GLUT_RIGHT_BUTTON:
        if(state == GLUT_DOWN){
            modifier = glutGetModifiers();
            previousMouse= vec2(x, y);
            if(modifier & GLUT_ACTIVE_CTRL){
                drag = ROTATE_CAMERA;
            }else{
                drag = SCALE;
            }
        }
        else if(state == GLUT_UP){
            drag = NO_DRAG;
        }
    break;
    default: break;
    }
    Context::display3DView();
}

void View3D::resetCamera(){
    Scene::eyeXRot = 0.f;
    Scene::eyeYRot = 0.f;

    Scene::eye = vec3(0,0,-8);
    updateCamera();
}

void View3D::updateCamera(){
    Scene::cameraMatrix = glm::lookAt(Scene::eye, vec3(0), vec3(0.0, 1.0, 0.0));
}
