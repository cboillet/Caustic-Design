/* ----------------------------------------------------------------
   name:           window.h
   purpose:        GLUT (sub-) window class declaration
   version:	   LIBRARY CODE
   author:         katrin lang
		   computer graphics
		   tu berlin
   ------------------------------------------------------------- */

#pragma once

#include <string>

class Window{

 public:

  // constructors
  Window();
  Window(const Window *parent, const std::string& title, int x, int y, int width, int height);

  void redisplay(void);

 protected:

  // make id current
  void select(void);

  // parent window
  Window *parent;

  // window identifier
  int id;

 public:

  // reshape window
  void reshape(int x, int y, int width, int height);

  //register callbacks with GLUT

  // register scene display function
  void registerDisplay(void display(void));
  
  // register scene redisplay function 
  // (after window reshape)
  void registerReshape(void reshape(int width, int height));
  
  // register mouse callback
  void registerMousePressed(void mousePressed(int btn, int state, int x, int y));
 
  // register mouse motion callback
  void registerMouseMoved(void mouseMoved(int x, int y));
 
  // register mouse motion callback
  void registerMouseDragged(void mouseMoved(int x, int y));

  // register keyboard callback
  void registerKeyPressed(void keyPressed(unsigned char key, int x, int y));
 
  // register keyboard callback for special keys
  void registerSpecialKeys(void specialKeys(int key, int x, int y));

  // register mouse menu
  void registerMenu(void menu(int id));
  
  // add mouse menu
  void addMenu(const char ids[], const std::string entries[], const int size);

  // add mouse menu
  void addMenu(const int ids[], const std::string entries[], const int size);

  // register idle callback (for animations)
  void registerIdle(void idle(void));    
};
