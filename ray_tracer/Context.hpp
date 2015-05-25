#ifndef CONTEXT_H
#define CONTEXT_H

#pragma once

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

   namespace Context{
    void init(int argc, char **argv);
    void printUsage();

    void display3DView();
    void displayRenderResult();
   }
#endif // CONTEXT_H
