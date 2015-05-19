/* ----------------------------------------------------------------
   name:           GLSLShader.cpp
   purpose:        cg1_ex3 OpenGL shading language wrapper
   version:	   LIBRARY CODE
   TODO:           nothing
   author:         katrin lang
                   computer graphics
                   tu berlin
   ------------------------------------------------------------- */

#include "GLSLShader.hpp"

using namespace std;

GLSLShader::GLSLShader(){

  program= 0;
}



GLSLShader::~GLSLShader(){

  for(unsigned int i= 0; i<vertexShaders.size(); i++)
    glDeleteShader(vertexShaders[i]);
  for(unsigned int i= 0; i<geometryShaders.size(); i++)
    glDeleteShader(geometryShaders[i]);
  for(unsigned int i= 0; i<fragmentShaders.size(); i++)
    glDeleteShader(fragmentShaders[i]);
  glDeleteProgram(program);
}

  // load vertex shader from file
void GLSLShader::loadVertexShader(const string& fileName){

  cout << "loading " << fileName << "...";
  compileVertexShader(load(fileName));
  if(vertexShaders.back()) cout << " done." << endl;
}
  // load geometry shader from file
void GLSLShader::loadGeometryShader(const string& fileName){

  cout << "loading " << fileName << "...";
  compileGeometryShader(load(fileName));
  if(geometryShaders.back()) cout << " done." << endl;
}
  // load fragment shader from file
void GLSLShader::loadFragmentShader(const string& fileName){

  cout << "loading " << fileName << "...";
  compileFragmentShader(load(fileName));
  if(fragmentShaders.back()) cout << " done." << endl;
}

  //  vertex shader directly from source contained in a string
void GLSLShader::compileVertexShader(const string& code){

  if(!program) program = glCreateProgram();
  GLuint id= compile(GL_VERTEX_SHADER, code);
  vertexShaders.push_back(id);
  if(id) glAttachShader(program, id);
}

// load geometry shader from directly from source contained in a string
void GLSLShader::compileGeometryShader(const string& code){

  if(!program) program = glCreateProgram();
  GLuint id= compile(GL_GEOMETRY_SHADER, code);
  geometryShaders.push_back(id);
  if(id) glAttachShader(program, id);
}
  // load fragment shader from directly from source contained in a string
void GLSLShader::compileFragmentShader(const string& code){

  if(!program) program = glCreateProgram();
  GLuint id= compile(GL_FRAGMENT_SHADER, code);
  fragmentShaders.push_back(id);
  if(id) glAttachShader(program, id);
}

void GLSLShader::link(void){
  //  cout << "linking... ";
  glLinkProgram(program);
  //cout << "done." << endl;
  printProgramLog(program);
  //TODO: invalidate hash tables
}

void GLSLShader::bind()const{

  glUseProgram(program);
}

void GLSLShader::unbind()const{

  glUseProgram(0);
}

void GLSLShader::setUniform(const string& name, const float& value){

  int id= getUniformBinding(name);
  glUniform1f(id, value);
}

void GLSLShader::setUniform(const string& name, const int& value){

  int id= getUniformBinding(name);
  glUniform1i(id, value);
}

void GLSLShader::setUniform(const string& name, const unsigned int& value){

  int id= getUniformBinding(name);
  glUniform1ui(id, value);
}

// Set a 2 component float vector uniform parameter
void GLSLShader::setUniform(const string& name, const float& x, const float& y){

  int id= getUniformBinding(name);
  glUniform2f(id, x, y);
}

// Set a 2 component integer vector uniform parameter
void GLSLShader::setUniform(const string& name, const int& x, const int& y){

  int id= getUniformBinding(name);
  glUniform2i(id, x, y);
}

// Set a 2 component undigned integer vector uniform parameter
void GLSLShader::setUniform(const string& name, const unsigned int& x, const unsigned int& y){

  int id= getUniformBinding(name);
  glUniform2ui(id, x, y);
}


// Set a 3 component float vector uniform parameter
void GLSLShader::setUniform(const string& name, const float& x, const float& y, const float& z){

  int id= getUniformBinding(name);
  glUniform3f(id, x, y, z);
}

// Set a 3 component integer vector uniform parameter
void GLSLShader::setUniform(const string& name, const int& x, const int& y, const int& z){

  int id= getUniformBinding(name);
  glUniform3i(id, x, y, z);
}

// Set a 3 component unsigned integer vector uniform parameter
void GLSLShader::setUniform(const string& name, const unsigned int& x, const unsigned int& y, const unsigned int& z){

  int id= getUniformBinding(name);
  glUniform3ui(id, x, y, z);
}

// Set a 4 component float vector uniform parameter
void GLSLShader::setUniform(const string& name, const float& x, const float& y, const float& z, const float& w){

  int id= getUniformBinding(name);
  glUniform4f(id, x, y, z, w);
}

// Set a 4 component integer vector uniform parameter
void GLSLShader::setUniform(const string& name, const int& x, const int& y, const int& z, const int& w){

  int id= getUniformBinding(name);
  glUniform4i(id, x, y, z, w);
}

// Set a 4 component unsigned integer vector uniform parameter
void GLSLShader::setUniform(const string& name, const unsigned int& x, const unsigned int& y, const unsigned int& z, const unsigned int& w){

  int id= getUniformBinding(name);
  glUniform4ui(id, x, y, z, w);
}

void GLSLShader::setUniform(const string& name, const glm::vec2& value){

  int id= getUniformBinding(name);
  glUniform2fv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::ivec2& value){

  int id= getUniformBinding(name);
  glUniform2iv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::uvec2& value){

  int id= getUniformBinding(name);
  glUniform2uiv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::vec3& value){

  int id= getUniformBinding(name);
  glUniform3fv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::ivec3& value){

  int id= getUniformBinding(name);
  glUniform3iv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::uvec3& value){

  int id= getUniformBinding(name);
  glUniform3uiv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::vec4& value){

  int id= getUniformBinding(name);
  glUniform4fv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::ivec4& value){

  int id= getUniformBinding(name);
  glUniform4iv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::uvec4& value){

  int id= getUniformBinding(name);
  glUniform4uiv(id, 1, &value[0]);
}

void GLSLShader::setUniform(const string& name, const glm::mat2& value){

  int id= getUniformBinding(name);
  glUniformMatrix2fv(id, 1, false, &value[0][0]);
}

void GLSLShader::setUniform(const string& name, const glm::mat3& value){

  int id= getUniformBinding(name);
  glUniformMatrix3fv(id, 1, false, &value[0][0]);
}

void GLSLShader::setUniform(const string& name, const glm::mat4& value){

  int id= getUniformBinding(name);
  glUniformMatrix4fv(id, 1, false, &value[0][0]);
}

void GLSLShader::bindVertexAttrib(const string& name, int id){

  glBindAttribLocation(program, id, name.c_str());
}

GLint GLSLShader::getVertexAttribBinding(const string& name){

  return glGetAttribLocation(program, name.c_str());
}

void  GLSLShader::bindOutput(const string& name, int id){
  glBindFragDataLocation(program, id, name.c_str());
}

GLint GLSLShader::getUniformBinding(const string& name){
  return glGetUniformLocation(program, name.c_str());
}

const string GLSLShader::load(const std::string& file){

  ifstream input(file.c_str(), ios::binary);
  assert(input.is_open());

  string line, code;

  while(getline(input,line)){
    code += line;
    code += "\n";
  }
  assert(!code.empty());
  return code;
}

GLuint GLSLShader::compile(GLuint type, const std::string& code){

  GLuint shader= glCreateShader(type);
  const char* codePointer= code.c_str();
  glShaderSource(shader, 1, &codePointer, 0);
  glCompileShader(shader);

  if(hasErrors(shader)){
    printShaderLog(shader);
    glDeleteShader(shader);
    return 0;
  }
  return shader;
}

bool GLSLShader::hasErrors(GLuint shader){
  GLint status;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  return status != GL_TRUE;
}

void GLSLShader::printProgramLog(GLuint program){

  GLint infoLogLength = 0;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);

  // Some drivers return 1 as infoLog_length when the infoLog is an empty string
  if(infoLogLength > 1) {
    char* infoLog = new char[infoLogLength];
    glGetProgramInfoLog(program, infoLogLength, 0, infoLog);
    cerr << "\nprogram log:" << endl << infoLog << endl << endl;
    delete [] infoLog;
  }
}

void GLSLShader::printShaderLog(GLuint shader){

  GLint infoLogLength = 0;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);

  // Some drivers return 1 as infoLog_length when the infoLog is an empty string
  if(infoLogLength > 1){
    char* infoLog = new char[infoLogLength];
    glGetShaderInfoLog(shader, infoLogLength, 0, infoLog);
    cerr << "\nshader log:" << endl << infoLog << endl << endl;
    delete [] infoLog;
  }
}
