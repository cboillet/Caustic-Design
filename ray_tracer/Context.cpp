#include "Context.hpp"
#include "Window.hpp"
#include "View3D.hpp"
#include "ViewRender.hpp"
#include "Scene.hpp"
#include "Assets.hpp"

#include <glm/glm.hpp>
#include <iostream>

#define GLM_FORCE_RADIANS

using namespace glm;
using namespace std;

// screen size
static uvec2 screen= uvec2(1024, 512);

// initial window position
static const uvec2 position= uvec2(100, 100);

// gap between subwindows
static const int GAP= 7;

//windows
static Window parentWindow, view3D, viewRender;



// display callback for GLUT
static void display(void){

  // clear color and depth buffer
  glClearColor(0.8, 0.8, 0.8, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glutSwapBuffers();
}

// reshape-callback for GLUT
static void reshape(int width, int height){

  // select main window
  glViewport(0, 0, width, height);

  screen= vec2(width, height);

  width-=3*GAP; height-=2*GAP;
  width/=2;

  // reshape 3d window
  view3D.reshape(GAP, GAP, width, height);

  // reshape render window
  viewRender.reshape(width+2*GAP, GAP, width, height);
}

static void specialKeys(int key, int x, int y){
    //unused variables
    (void) x;
    (void) y;

    switch(key){
    case GLUT_KEY_F12:
        Scene::render();
        break;

    case GLUT_KEY_F10:
        Scene::texturingEnabled = (!Scene::texturingEnabled);
        cout << "Texturing enabled: ";
        if(Scene::texturingEnabled) cout << "true" << endl;
        else cout << "false" << endl;
        break;

    case GLUT_KEY_UP:
        Scene::recursionDepth++;
        cout << "Recursion-Depth = " << Scene::recursionDepth << endl;
        break;

    case GLUT_KEY_DOWN:
        Scene::recursionDepth--;
        if(Scene::recursionDepth < 0){
            Scene::recursionDepth = 0;
        }
        cout << "Recursion-Depth = " << Scene::recursionDepth << endl;
     default:
        break;

    }
}



//creating main window, subwindow plus registering callbacks
static void createWindows(){
    parentWindow = Window(NULL, "Raytracer", position.x, position.y, screen.x, screen.y);
    parentWindow.registerDisplay(display);
    parentWindow.registerReshape(reshape);
    parentWindow.registerSpecialKeys(specialKeys);

    int subWidth= (screen.x - 3*GAP)/2;
    int subHeight= screen.y - 2*GAP;

    view3D= Window(&parentWindow, "", GAP, GAP, subWidth, subHeight);
    view3D.registerDisplay(View3D::display);
    view3D.registerReshape(View3D::reshape);
    view3D.registerKeyPressed(View3D::keyPressed);
    view3D.registerMouseDragged(View3D::mouseDragged);
    view3D.registerMousePressed(View3D::mousePressed);
    view3D.registerSpecialKeys(specialKeys);


    //TODO this is only for testing. should be through ui later
    Scene::addNewModel();
    Scene::setLighting();
    //Scene::setMaterial();

    viewRender= Window(&parentWindow, "", subWidth + 2*GAP, GAP, subWidth, subHeight);
    viewRender.registerDisplay(ViewRender::display);
    viewRender.registerReshape(ViewRender::reshape);
    viewRender.registerKeyPressed(ViewRender::keyPressed);
    viewRender.registerMouseDragged(ViewRender::mouseDragged);
    viewRender.registerMousePressed(ViewRender::mousePressed);
    viewRender.registerSpecialKeys(specialKeys);

}

void Context::display3DView() {
    view3D.redisplay();
}

void Context::displayRenderResult() {
    viewRender.redisplay();
}

// initialize OpenGL context
void Context::init(int argc, char **argv){

  // create window with glut
  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

  glutInitContextVersion(3,3);
  glutInitContextProfile(GLUT_CORE_PROFILE);


  glutSetOption(GLUT_RENDERING_CONTEXT, GLUT_USE_CURRENT_CONTEXT);


  printUsage();

  createWindows();

  glewExperimental = true;
GLenum err= glewInit();
#ifndef __APPLE__
  if(err != GLEW_OK){
     cerr << "GLEW not available!" << endl;
  }
   #endif

  if(GLEW_VERSION_3_3)
      cout << "version 3.3 supported"<< endl;

  glGetError();


  View3D::loadShader();
  ViewRender::init();

  cout << "GPU: " << glGetString(GL_RENDERER) << ", OpenGL version: " << glGetString(GL_VERSION) << endl;


 glutMainLoop();
}

void Context::printUsage(){
    cout << "=================================================================" << endl;
    cout << "USAGE:" << endl;
    cout << "----------GENERAL------------------------------------------------" << endl;
    cout << "q/Q:\t\t Exit" << endl;
    cout << "F12:\t\t Start Ray Tracer" << endl << endl;
    cout << "----------LEFT WINDOW--------------------------------------------" << endl;
    cout << "Left-Mouse:\t\tRotate selected object" << endl;
    cout << "Ctrl+Left-Mouse:\tTranslate object in xy-plane" << endl;
    cout << "Shift+Left-Mouse:\tTranslate object in z-direction" << endl;
    cout << "Right-Mouse:\tScale Object" << endl;
    cout << "Ctrl+Right-Mouse:\tRotate Camera around Origin" << endl;
    cout << "m/M:\t\tSkip through models" << endl;
    cout << "e/E:\t\tSkip through elements (objects)" << endl;
    cout << "t/T:\t\t+/- red-channel of light-source" << endl;
    cout << "g/G:\t\t+/- green-channel of light-source" << endl;
    cout << "b/B:\t\t+/- blue-channel of light-source" << endl;
    cout << "r/R:\t\tReset Camera Position" << endl;
    cout << "n: \t\tAdd new element" << endl;
    cout << "d: \t\tDelete current element" << endl << endl;
    cout << "----------RIGHT WINDOW-------------------------------------------" << endl;
    cout << "Left-Mouse:\t\tRotate Camera Around Origin" << endl;
    cout << "r/R:\t\tReset Camera Position" << endl << endl;
    cout << "=================================================================" << endl;
}





int main(int argc, char** argv){

  Context::init(argc, argv);

  return 0;
}
