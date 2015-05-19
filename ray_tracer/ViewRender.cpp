#include "ViewRender.hpp"
#include "Scene.hpp"
#include "Context.hpp"
#include "TriMesh.hpp"
#include "GLSLShader.hpp"

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

#define PI  3.14159265358979323846264338327950288f
#define RADIANS(x) (((x)*PI)/180.0f)

using namespace glm;
using namespace std;

TriMesh quad("meshes/quad.off", false);

bool buffered = false;

// current state of mouse action
static enum{
    ROTATE, SHIFT_XY, SHIFT_Z, SCALE, NO_DRAG, DRAW, ERASE, SCALEBRUSH, SMOOTHBRUSH
} drag= NO_DRAG;

vec2 prevMouse = vec2();

std::vector<glm::vec3> points;
std::vector<glm::vec3> colors;

GLSLShader imageViewShader;

vec2 ViewRender::screen;

float eyeDistance = -3.0f;
float eyeXRot = 0.0f;
float eyeYRot = 0.0f;
vec3 initialEyePos = vec3(0.0, 0.0, -8.0);
vec3 eye = vec3(initialEyePos);
mat4 cameraMatrix = glm::lookAt(eye, vec3(0), vec3(0.0, 1.0, 0.0));


GLuint vao = 0;

void ViewRender::init(){
    imageViewShader.loadVertexShader("shaders/view3d.vert");
    imageViewShader.loadFragmentShader("shaders/view3d.frag");
    imageViewShader.bindVertexAttrib("position", TriMesh::attribVertex);
    imageViewShader.bindVertexAttrib("color", TriMesh::attribColor);
    //imageViewShader.bindVertexAttrib("texcoord", TriMesh::attribTexCoord);
    imageViewShader.link();

    Scene::image.setBlack(Scene::xRes, Scene::yRes);
    Scene::image.generateTexture();
    quad.generateDefaultUV();
}

void ViewRender::buffer(){

    if(vao != 0){
        glDeleteVertexArrays(1, &vao);
    }

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLuint vbos[2];
    glGenBuffers(2, &vbos[0]);

    glBindBuffer(GL_ARRAY_BUFFER, vbos[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*points.size(), &points[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(TriMesh::attribVertex);
    glVertexAttribPointer(TriMesh::attribVertex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*colors.size(), &colors[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(TriMesh::attribColor);
    glVertexAttribPointer(TriMesh::attribColor, 3, GL_FLOAT, GL_FALSE, 0, 0);

    buffered = true;

    glBindVertexArray(0);
}

void ViewRender::display(){

    if(!buffered){
        //normalizeColors();
        buffer();
    }

    glClearColor(0.5, 0.5, 0.5, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    imageViewShader.bind();
    imageViewShader.setUniform("viewProjectionMatrix", Scene::projectionMatrix*cameraMatrix);

    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, points.size());
    glBindVertexArray(0);
    imageViewShader.unbind();


    glutSwapBuffers();

}

void ViewRender::reshape(int width, int height){
    glViewport(0, 0, width, height);

    screen = vec2(width, height);
}

void ViewRender::mouseDragged(int x, int y){

    glm::mat4 rotX;
    glm::mat4 rotY;
    vec2 v= (vec2(x,y) - prevMouse) / screen;
    switch(drag){
        case ROTATE:
            eyeYRot -= v.x*3.f;
            eyeXRot += v.y*3.f;
            rotX = glm::rotate(mat4(1), eyeYRot, vec3(0,1,0));
            rotY = glm::rotate(mat4(1), eyeXRot, vec3(1,0,0));
            eye = vec3(rotX * rotY * vec4(initialEyePos,1));
            updateCamera();
            break;
        case SHIFT_XY:
            //Scene::translate (-vec3(3.3*v.x, 3.3*v.y, 0));
            break;
        case SHIFT_Z:
            //Scene::translate(vec3(0.0, 0.0, 3.3*sign(dot(v, vec2(1,1))) * length(v)));
            break;
        default:
            break;
    }

    prevMouse= vec2(x, y);
    Context::displayRenderResult();
}

void ViewRender::mousePressed(int button, int state, int x, int y){

    int modifier;

    switch(button){
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN){
        prevMouse= vec2(x, y);
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
    default: break;
    }
    Context::displayRenderResult();
}

void ViewRender::keyPressed(unsigned char key, int x, int y){

    switch (key){
        case 'q':
        case 'Q':
            exit(EXIT_SUCCESS);

        case 'r':
            resetCamera();
            break;
        case 'R':
            resetCamera();
            break;
    }

    Context::displayRenderResult();
}

void ViewRender::updateImage(){
    Scene::image.generateTexture();
}

void ViewRender::addPoint(vec3 pointPosition, vec3 pointColor){
    //float x = pointPosition.x, y = pointPosition.y, z = pointPosition.z;
    //cout << "adding point at (" << pointPosition.x << ", " << pointPosition.y << ", " << pointPosition.z << ")" << endl;
    //cout << "color is (" << pointColor.x << ", " << pointColor.y << ", " << pointColor.z << ")" << endl;
    points.push_back(pointPosition);
    colors.push_back(pointColor);

    buffered = false;
}

void ViewRender::normalize(){

    float max = 0.f;
    for(int i=0; i<colors.size(); i++){
        for(int j=0; j<3; j++){
            if(colors[i][j] > max){
                max = colors[i][j];
            }
        }
    }

    if(max > 1.0f){
        for(int i=0; i<colors.size(); i++){
            for(int j=0; j<3; j++){
                colors[i][j] /= max;
            }
        }
    }
    buffered = false;
}

void ViewRender::removeOldPoints(){
    points.clear();
    colors.clear();

    buffered = false;
}

void ViewRender::resetCamera(){
    eyeXRot = Scene::eyeXRot;
    eyeYRot = Scene::eyeYRot;

    eye = Scene::eye;
    updateCamera();
}

void ViewRender::updateCamera(){
    cameraMatrix = glm::lookAt(eye, vec3(0), vec3(0.0, 1.0, 0.0));
}
