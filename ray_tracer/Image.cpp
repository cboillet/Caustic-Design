/* ----------------------------------------------------------------
   name:           Image.cpp
   purpose:        cg1_ex4 ws2014/15 texturing tutorial
   version:	   SKELETON CODE
   TODO:           texture and mipmap generation, texture filtering, wrapping, texel get, painting in texture (see XXX)
   author:         katrin lang
                   computer graphics
                   tu berlin
   ------------------------------------------------------------- */

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#ifdef WITH_OPENIMAGEIO
	#include <OpenImageIO/imageio.h>
#endif

#include "Image.hpp"
#include "Context.hpp"

using namespace std;
using namespace glm;

#ifdef WITH_OPENIMAGEIO
	OIIO_NAMESPACE_USING
#endif

Image::Image() : width(0), height(0), textureID(0), wrapS(GL_REPEAT), wrapT(GL_REPEAT), mag(GL_LINEAR), min(GL_LINEAR),  trilinear(false){
}

Image::Image(int width, int height)
  : data(width*height)
  , width(width)
  , height(height)
  , textureID(0)
  , wrapS(GL_REPEAT)
  , wrapT(GL_REPEAT)
  , mag(GL_LINEAR)
  , min(GL_LINEAR)
  , trilinear(false)
{
}

Image::Image(const std::string& filename) : textureID(0), wrapS(GL_REPEAT), wrapT(GL_REPEAT), mag(GL_LINEAR), min(GL_LINEAR), trilinear(false){
  load(filename);
}

Image::~Image(){
}

void Image::normalize(){
    float max = 0.f;

    for(int i=0; i<data.size(); i++){
        for(int j=0; j<3; j++){
            if(data[i][j] > max){
                max = data[i][j];
            }
        }
    }

    if(max > 1.0f){
        for(int i=0; i<data.size(); i++){
            for(int j=0; j<3; j++){
                data[i][j] /= max;
            }
        }
    }

}

void Image::setBlack(int width, int height){
    data.resize(width*height);

    for(int i=0; i<width*height; i++){
        data[i] = vec4(0, 0, 0, 1);
    }
}

// generate OpenGL texture
// XXX: NEEDS TO BE IMPLEMENTED
void Image::generateTexture(){

  if(textureID != 0){
      glDeleteTextures(1, &textureID);
      textureID = 0;
  }

  if(textureID==0){
    // generate texture id
    // XXX
      glGenTextures(1, &textureID);
    // END XXX
  }

  bind();
  // texture filtering and repeat
  // XXX

  setMagFilter(mag);
  setMinFilter(min);

  setWrapS(wrapS);
  setWrapT(wrapT);
  // END XXX

  //enable automatic mipmap generation
  // XXX
  if(trilinear){
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_FALSE);
  }else{
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
  }
  // END XXX

  // upload texture data
  // XXX
  if(trilinear){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (int)((float)width)*1.5, height, 0, GL_RGBA, GL_FLOAT, &mipmap[0]);
  }else{
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, &data[0]);
  }
  //glGenerateMipmap(GL_TEXTURE_2D);
  // END XXX
  unbind();
}

void Image::setMinFilter(GLuint min){
  this->min= min;

  // set texture parameter
  // XXX

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min);

  // END XXX
}

// set magnifying filter
// XXX: NEEDS TO BE IMPLEMENTED
void Image::setMagFilter(GLuint mag){

  this->mag= mag;

  // set texture parameter
  // XXX

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag);

  // END XXX
}

// set wrapping mode
// XXX: NEEDS TO BE IMPLEMENTED
void Image::setWrapS(GLuint wrap){

  this->wrapS= wrap;

  // set texture filter
  // XXX
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapS);

  // END XXX
}

// set wrapping mode
// XXX: NEEDS TO BE IMPLEMENTED
void Image::setWrapT(GLuint wrap){

  this->wrapT= wrap;

  // set texture filter
  // XXX
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapT);

  // END XXX
}

// set both wrapping modes
void Image::setWrap(GLuint wrap){
  setWrapS(wrap);
  setWrapT(wrap);
}

void Image::setTrilinear(bool trilin){
    this->trilinear = trilin;
    generateTexture();
}


// bind texture
// XXX: NEEDS TO BE IMPLEMENTED
void Image::bind(){
  // bind texture
  // XXX
    glBindTexture(GL_TEXTURE_2D, textureID);
    //cout << "binding texture with id " << textureID << endl;
  // END XXX
}

// unbind texture
// XXX: NEEDS TO BE IMPLEMENTED
void Image::unbind(){
  // XXX
    glBindTexture(GL_TEXTURE_2D, 0);
  // END XXX
}

void Image::generateMipMap(){

    // NOTE:	using mipmap as given here: http://de.wikipedia.org/wiki/Datei:MipMap_Example_STS101.jpg

    // image in original size
    for(unsigned int i=0; i < data.size(); i++){
        int x = i % (width);
        int y = i / (width);

        mipmapdata[x][y] = data[i];
    }

    //if(false){
    // sub-parts

    // dimensions of current mipmap-level we write
    int currentWidth = width / 2;
    int currentHeight = height / 2;
    // defines origin of already existing parent image
    int readOffsetX = 0;
    int readOffsetY = 0;
    // origin of image we create
    int writeOffsetX = width;
    int writeOffsetY = 0;

    int lod = 1;
    while(currentWidth >= 1 && currentHeight >= 1){
        for(int u=0; u<currentWidth; u++){
            for(int v=0; v<currentHeight; v++){
                int readIndexX = readOffsetX + u*2;
                int readIndexY = readOffsetY + v*2;
                int writeIndexX = writeOffsetX + u;
                int writeIndexY = writeOffsetY + v;

                mipmapdata[writeIndexX][writeIndexY] = mipmapdata[readIndexX][readIndexY];
                mipmapdata[writeIndexX][writeIndexY] += mipmapdata[readIndexX+1][readIndexY];
                mipmapdata[writeIndexX][writeIndexY] += mipmapdata[readIndexX][readIndexY+1];
                mipmapdata[writeIndexX][writeIndexY] += mipmapdata[readIndexX+1][readIndexY+1];
                mipmapdata[writeIndexX][writeIndexY] /= 4.0;

            }
        }

        if(lod == 1){
            readOffsetX = writeOffsetX;
        }

        readOffsetY = writeOffsetY;
        writeOffsetY += currentHeight;


        currentWidth /= 2;
        currentHeight /= 2;
        lod++;
    }
    //}
    for(unsigned int i=0; i<mipmap.size(); i++){
        int x = i % (int)(1.5*width);
        int y = i / (int)(1.5*width);
        mipmap[i] = mipmapdata[x][y];
    }

    // debugging output -- outcoming image should look like original (and it does)
    // for further testing, smaller mipmap-levels should be used
    /*
    xOffset = width;
    yOffset = height;
    for (int x = 0; x < width; x++){
        for(int y = 0; y < height; y++){
            data[x+y*width] = vec4(mipmap[x][y], mipmap[x+xOffset][y], mipmap[x][y+yOffset], 1.0);
        }
    }*/
}

void Image::getPatch(unsigned int xoffset, unsigned int yoffset, unsigned int patchWidth, unsigned int patchHeight, vector<vec4> &patch){
    if(patch.size() == patchWidth*patchHeight){
        int upperLeft = (yoffset*this->width)+xoffset;
        for(unsigned int i = 0; i < patchWidth*patchHeight; i++){
           patch[i] = this->data[upperLeft+ i%patchWidth + (i/patchWidth)*this->width];
        }

    }else{
        cerr << "ERROR: getPatch was called with wrong paramters!"<<endl;
    }
}

// get texture pixels
std::vector<glm::vec4>* Image::getPixels(){
    return &this->data;
}

unsigned int Image::getHeight(){
    return this->height;
}

unsigned int Image::getWidth(){
    return this->width;
}

void Image::setSize(unsigned int width, unsigned int height){
    data.clear();
    data.resize(width*height);
    this->width = width;
    this->height = height;
}

void Image::load(const std::string& filename){

  name = filename;

  data.clear();

  if(filename.substr(filename.size()-4, 4) == ".ppm"){
      loadPPM(filename);
      generateMipMap();
  }
  else{
    cerr << "file " << filename << " is not a PPM file" << endl;
    return;
  }
}

void Image::loadPPM(const std::string& filename){

  ifstream file(filename.c_str(), ios::binary);

  if(!file.is_open()){
    cerr << "opening file " << filename << " failed" << endl;
    return;
  }

  // grab first two chars of the file and make sure that it has the
  // correct magic cookie for a raw PPM file.
  string magic;
  getline(file, magic);
  if(magic.substr(0, 2) != "P6"){
    cerr << "File " << filename << " is not a raw PPM file" << endl;
    return;
  }

  // grab the three elements in the header (width, height, maxval).
  string dimensions;
  do{
    getline(file, dimensions);
  }
  while(dimensions[0] == '#');

  stringstream(dimensions) >> width >> height;

  string max;
  getline(file, max);
  int maxValue;
  stringstream(max) >> maxValue;
  // grab all the image data in one fell swoop.
  vector<char> raw(width*height*3);
  file.read(&raw[0], raw.capacity());
  file.close();

  data.resize(width*height);
  mipmap.resize(width*1.5*height);

  mipmapdata.resize(width*1.5);
  for(unsigned int i=0; i<mipmapdata.size(); i++){
      mipmapdata[i].resize(height);
  }


  for(int y = 0; y < height; y++){
    for(int x = 0; x < width; x++){
      data[y*width+x]= vec4((unsigned char)raw[(height - y-1) * width * 3 + 3*x], (unsigned char)raw[(height - y-1) * width * 3 + 3*x + 1], (unsigned char)raw[(height - y-1) * width * 3 + 3*x + 2], maxValue);
      data[y*width+x]/= maxValue;
      //cout << data[i].r << " " + data[i].g << " " + data[i].b << " " + data[i].a << endl;
    }
  }

  raw.clear();

  std::cout << "Image " << filename << " loaded. width=" << width << " height=" << height << endl;
}

void Image::setPixel(unsigned int  index, glm::vec3 color){
    this->data[index] = vec4(color, 1.0);
}

vec4 Image::getRelativePixel(vec2 &uv){
    return getPixel(int(uv.x*float(width)), int(uv.y*float(height)));
}

vec4 Image::getPixel(unsigned int x, unsigned int y){
    x = x % width;
    y = y % height;
    return this->data[this->width*y+x];
}

void Image::savePPM(){

    // normalize values
    float max = 0.0f;
    for(unsigned int i=data.size(); i--;){
        for(unsigned int j=3; j--;){
            if(data[i][j] > max){
                max = data[i][j];
            }
        }
    }
    for(unsigned int i=data.size(); i--;){
        for(unsigned int j=3; j--;){
            data[i][j] /= max;
        }
    }



    ofstream file;
    file.open("render.ppm", ios::out | ios::binary);
    if(file.is_open()){

        // magic numer
        file << "P6" << endl;
        file << width << " " << height << endl;
        file << "255" << endl;

        for(unsigned int y=height; y--;){
            for(unsigned int x = 0; x<width; x++){
                unsigned int i = x + y * width;
                for(unsigned int j=0; j<3; j++){
                    file << static_cast<char>((int)(255.0*data[i][j]));
                }
            }
        }

        file.close();
    }else{
        cerr << "ERROR: Could not open file render.ppm" << endl;
    }
}

#ifdef WITH_OPENIMAGEIO
void Image::saveToDisk(string name){
    ImageOutput *out = ImageOutput::create(name);
    if(!out){
        cerr << "could not write image"<< endl;
        return;
    }
    ImageSpec spec (width, height, 4, TypeDesc::FLOAT);

    out->open(name, spec);
    out->write_image(TypeDesc::FLOAT, this->data.data());
    out->close();

    delete out;
}
#endif


