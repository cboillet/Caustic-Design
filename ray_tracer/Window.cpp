/* ----------------------------------------------------------------
   name:           Window.cpp
   purpose:        GLUT (sub-) window implementation
   version:	   LIBRARY CODE
   author:         katrin lang
		   computer graphics
		   tu berlin
   ------------------------------------------------------------- */

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

#include "Window.hpp"

using namespace std;

Window::Window() : parent(NULL) {}

Window::Window(const Window *parent, const string& title, int x, int y, int w, int h){

  if(!parent){ // main window
    glutInitWindowPosition(x, y);
    glutInitWindowSize(w, h);
    id= glutCreateWindow(title.c_str());
  }
  else id= glutCreateSubWindow(parent->id, x, y, w, h);
}

void Window::select(void){
  glutSetWindow(id);
}

void Window::redisplay(void){
  select();
  glutPostRedisplay();
}

  // reshape window
void Window::reshape(int x, int y, int width, int height){

  select();
  glutPositionWindow(x, y);
  glutReshapeWindow(width, height);
  glutPostRedisplay();
}

//register callbacks with GLUT

  // register scene display function
void Window::registerDisplay(void display(void)){
  select();
  glutDisplayFunc(display);
}

// register scene redisplay function 
// (after window reshape)
void Window::registerReshape(void reshape(int width, int height)){
  select();
  glutReshapeFunc(reshape);
}

// register mouse callback
void Window::registerMousePressed(void mousePressed(int btn, int state, int x, int y)){
  select();
  glutMouseFunc(mousePressed);
}

// register mouse motion callback
void Window::registerMouseMoved(void mouseMoved(int x, int y)){
  select();
  glutPassiveMotionFunc(mouseMoved);
}

// register mouse motion callback
void Window::registerMouseDragged(void mouseMoved(int x, int y)){
  select();
  glutMotionFunc(mouseMoved);
}

// register keyboard callback
void Window::registerKeyPressed(void keyPressed(unsigned char key, int x, int y)){
  select();
  glutKeyboardFunc(keyPressed);
}

// register keyboard callback for special keys
void Window::registerSpecialKeys(void specialKeys(int key, int x, int y)){
  select();
  glutSpecialFunc(specialKeys);
}

// register mouse menu
void Window::registerMenu(void menu(int id)){
  select();
  glutCreateMenu(menu);
}

// fill in mouse menu
void Window::addMenu(const char ids[], const string entries[], const int n){

  select();
  for(int i= 0; i<n; i++){
    glutAddMenuEntry(entries[i].c_str(), ids[i]);
  }
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// fill in mouse menu
void Window::addMenu(const int ids[], const string entries[], const int n){

  select();
  for(int i= 0; i<n; i++){
    glutAddMenuEntry(entries[i].c_str(), ids[i]);
  }
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

  // register idle callback (for animations)
void Window::registerIdle(void idle(void)){
  select();
  glutIdleFunc(idle);
}
