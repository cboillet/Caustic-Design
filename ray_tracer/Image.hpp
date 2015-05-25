/* ----------------------------------------------------------------
   name:           Image.hpp
   purpose:        cg1_ex4 ws2014/15 texturing tutorial
   version:	   SKELETON CODE
   TODO:           nothing (see Image.cpp)
   author:         katrin lang
                   computer graphics
                   tu berlin
   ------------------------------------------------------------- */

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

#include <vector>
#include <string>

#include "glm/glm.hpp"

class Image{

public:

  // constructors
  Image();
  Image(int width, int height);
  Image(const std::string& filename);

  // destructor
  ~Image();

  void setBlack(int width, int height);

  std::string getName(){ return name;}
  // load image from file
  void load(const std::string& filename);

  // set texture filter
  void setMinFilter(GLuint min);
  void setMagFilter(GLuint mag);

  // set wrapping mode
  void setWrapS(GLuint wrap);
  void setWrapT(GLuint wrap);
  // set both S and T
  void setWrap(GLuint wrap);

  void setTrilinear(bool trilin);

  // bind/unbind texture
  void bind();
  void unbind();

  // generate OpenGL texture
  void generateTexture();

  void normalize();

  // draw in texture
  void setPixel(unsigned int  index, glm::vec3 color);

  unsigned int getHeight();
  unsigned int getWidth();
  void setSize(unsigned int width, unsigned int height);

  void getPatch(unsigned int xoffset, unsigned int yoffset, unsigned int patchWidth, unsigned int patchHeight, std::vector<glm::vec4> &patch);

  std::vector<glm::vec4>* getPixels();
  glm::vec4 getPixel(unsigned int x, unsigned int y);
  glm::vec4 getRelativePixel(glm::vec2 & uv);

  void savePPM();
  void saveToDisk(std::string name);


protected:

  // mipmap
  //std::vector<glm::vec4> mipmap;
  std::vector<glm::vec4 > mipmap;
  std::vector<std::vector<glm::vec4> > mipmapdata;

  // image data
  std::vector<glm::vec4> data;
  // dimensions
  int width;
  int height;

  std::string name;

  GLuint textureID;

  //texturing parameters
  GLuint wrapS;
  GLuint wrapT;
  GLuint mag;
  GLuint min;

  bool trilinear;

  int brushSize;
  std::vector<glm::vec4> brush;

  // parse ppm format
  void loadPPM(const std::string& filename);

  void generateMipMap(void);

};

