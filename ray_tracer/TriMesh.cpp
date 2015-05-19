/* ----------------------------------------------------------------
   name:           TriMesh.cpp
   purpose:        cg1_ex3 2014 triangle mesh for OpenGL rendering
   version:	   SKELETON CODE
   author:         katrin lang
                   computer graphics
                   tu berlin
   ------------------------------------------------------------- */

#include "TriMesh.hpp"

// use this with care
// might cause name collisions
using namespace glm;

using namespace std;

// NVIDIA wants it like this
// see https://www.opengl.org/sdk/docs/tutorials/ClockworkCoders/attributes.php
const GLuint TriMesh::attribVertex= 0;
const GLuint TriMesh::attribNormal= 2;
const GLuint TriMesh::attribColor= 3;
const GLuint TriMesh::attribTexCoord= 8;


#define PI  3.14159265358979323846264338327950288f
#define RADIANS(x) (((x)*PI)/180.0f)



TriMesh::TriMesh(){
  winding= CW;
    vaos = vector<GLuint>(2);
    selected = false;
}

TriMesh::TriMesh(const std::string& fileName){
    vaos = vector<GLuint>(2);
    selected = false;
name= fileName;
  winding= CW;
  loadOff(fileName);
  center();
  unitize();
  computeNormals();
}


TriMesh::TriMesh(const std::string& fileName, bool doUnitize){
    vaos = vector<GLuint>(2);
    selected = false;
    name = fileName;
    winding = CW;
    cout << "before load"<< endl;
    loadOff(fileName);
    cout << "loaded"<< endl;
    if(doUnitize){
        center();
        unitize();
    }
    computeNormals();
}

TriMesh::~TriMesh(){
    vaos = vector<GLuint>(2);
    selected = false;
}

void TriMesh::setWinding(PolygonWinding winding){
  this->winding= winding;
}

// center model at its origin
void TriMesh::center(void){

  calculateBoundingBox();

  vec3 center= (boundingBoxMin + boundingBoxMax) * vec3(0.5);

  for(unsigned int i= 0; i<positions.size(); i++){
    positions[i]-= center;
  }
  boundingBoxMin-= center;
  boundingBoxMax-= center;
}

// normalize to bounding sphere radius 1
void TriMesh::unitize(void){

  calculateBoundingSphere();

  for(unsigned int i= 0; i<positions.size(); ++i){
    positions[i]/= boundingSphereRadius;
  }
  boundingSphereRadius= 1;
  boundingBoxMin= vec3(-1);
  boundingBoxMax= vec3(1);
}

// calculate bounding sphere
void TriMesh::calculateBoundingSphere(void){

  boundingSphereRadius= 0;
  for(unsigned int i= 0; i<positions.size(); i++){
    vec3 v= positions[i];
    if(length(v) > boundingSphereRadius) boundingSphereRadius= length(v);
  }
}

// calculate bounding box
void TriMesh::calculateBoundingBox(void){

  boundingBoxMin= vec3(numeric_limits<float>::max());
  boundingBoxMax= vec3(numeric_limits<float>::min());
  for(unsigned int i= 0; i<positions.size(); i++){
    if(positions[i].x < boundingBoxMin.x) boundingBoxMin.x= positions[i].x;
    if(positions[i].x > boundingBoxMax.x) boundingBoxMax.x= positions[i].x;
    if(positions[i].y < boundingBoxMin.y) boundingBoxMin.y= positions[i].y;
    if(positions[i].y > boundingBoxMax.y) boundingBoxMax.y= positions[i].y;
    if(positions[i].z < boundingBoxMin.z) boundingBoxMin.z= positions[i].z;
    if(positions[i].z > boundingBoxMax.z) boundingBoxMax.z= positions[i].z;
  }
}

void TriMesh::correctTexture(bool correct){
  textureCorrection= correct;
}

// load triangle mesh in OFF format
void TriMesh::reload(){
  loadOff(name);
}

// load triangle mesh in .OFF format
void TriMesh::loadOff(const string& fileName){

    deleteVAO();

    name= fileName;

    positions.clear();
    faces.clear();
    normals.clear();

    ifstream file;
    file.open(fileName.c_str());

    string tmp;

    getline(file, tmp);
    /*if(line.compare("OFF") != 0){
        //Invalid header
    }*/

    getline(file, tmp);

    stringstream line(tmp);
    unsigned int v;
    unsigned int f;
    unsigned int e;

    line >> v >> f >> e;

    for(unsigned int i = 0; i < v; i++){
        getline(file, tmp);
        stringstream line(tmp);

        float x, y, z;

        line >> x >> y >> z;

        positions.push_back(vec3(x, y, z));

        //cout << "found position: "<<x << ", " <<y  << ", "  <<z<<endl;
    }
    //cout<< "\n==============================\n"<<endl;

    for(unsigned int i = 0; i < f; i++){
        getline(file, tmp);
        stringstream line(tmp);

        int num, i1, i2, i3;

        line >> num >> i1 >> i2 >> i3;

        if(winding == CW){
            faces.push_back(uvec3(i1, i2, i3));
        }else{
            faces.push_back(uvec3(i3, i2, i1));
        }

    }


}

void TriMesh::getFaceNormal(unsigned int faceIndex, vec3 &normal){
    vec3 v1 = positions[faces[faceIndex][1]] - positions[faces[faceIndex][0]];
    vec3 v2 = positions[faces[faceIndex][2]] - positions[faces[faceIndex][0]];

    normal = normalize(cross(v1, v2));
}

// calculate smooth per-vertex normals
void TriMesh::computeNormals(void){
    vector<vec3> faceNormals;

    faceNormals.clear();
    //calculating normals per face
    for(unsigned int i = 0; i<faces.size(); i++){
        vec3 normal;
        getFaceNormal(i,normal);
        fflush(stdout);
        faceNormals.push_back(normal);
    }


    //aF contains indices of adjacend faces per vertex
    //vector<unsigned int>aF[positions.size()];
    vector<vector<unsigned int> > aF;
    aF.resize(positions.size());
    for(unsigned i = 0; i < faces.size(); i++){
        aF[faces[i][0]].push_back(i);
        aF[faces[i][1]].push_back(i);
        aF[faces[i][2]].push_back(i);
    }


    //interpolate normals of adjacend faces per vertex
    for(unsigned int i = 0; i< aF.size(); i++){
        vec3 vertexNormal = vec3(0);
        for(unsigned int j= 0; j < aF[i].size(); j++){
                vertexNormal += faceNormals[aF[i][j]];
                //cout << "current normal[" << j << "] is " << vertexNormal[0] << ", " << vertexNormal[1] << ", " << vertexNormal[2] << endl;
        }
        vertexNormal /= aF[i].size();
        //cout << "final normal is " << vertexNormal[0] << ", " << vertexNormal[1] << ", " << vertexNormal[2] << endl;

        normals.push_back(vertexNormal);
    }
}

  // Compute uv coordinates with a spherical mapping
  // (vertices are projected on a sphere along the normal and classical sphere uv unwrap is used)
void TriMesh::computeSphereUVs(void){
    texCoords.clear();

    float u, v; // theta is between y-axis and vector, phi is between x-axis and vector
    vec2 vec;
    vec3 upVec;
    //float minx = 1000000.0, miny = 10000000.0, maxx = -1.0, maxy = -1.0;
    for (unsigned int i=0; i<positions.size(); i++){
        upVec = normalize(vec3(positions[i])) * boundingSphereRadius;
        v = acos( -upVec.y );
        v /= PI;
        u = atan2(upVec.x, upVec.z);
        u /= 2*PI;
        u += 0.5;
        vec = vec2(u, v);

        texCoords.push_back(vec);
    }

    // do checking for texture problems
    cout << "correcting texture" << endl;
    for(unsigned int i=0; i<faces.size(); i++){

        vec2 uv1 = texCoords[faces[i][0]];
        vec2 uv2 = texCoords[faces[i][1]];
        vec2 uv3 = texCoords[faces[i][2]];

        // check u-coordinate
        bool small1 = uv1[0] < 0.25;
        bool small2 = uv2[0] < 0.25;
        bool small3 = uv3[0] < 0.25;

        bool big1 = uv1[0] > 0.75;
        bool big2 = uv2[0] > 0.75;
        bool big3 = uv3[0] > 0.75;

        vec2 texExtra(-1, -1);
        int index = -1;
        if(small1 && small2 && big3){
            index = 2;
            texExtra = vec2(uv3[0]-1.0, uv3[1]);
        }else if(small1 && big2 && small3){
            index = 1;
            texExtra = vec2(uv2[0]-1.0, uv2[1]);
        }else if(big1 && small2 && small3){
            index = 0;
            texExtra = vec2(uv1[0]-1.0, uv1[1]);
        }else if(big1 && big2 && small3){
            index = 2;
            texExtra = vec2(uv3[0]+1.0, uv3[1]);
        }else if(big1 && small2 && big3){
            index = 1;
            texExtra = vec2(uv2[0]+1.0, uv2[1]);
        }else if(small1 && big2 && big3){
            index = 0;
            texExtra = vec2(uv1[0]+1.0, uv1[1]);
        }

        if(index != -1){
            // new position is same as old position of given vertex
            vec3 v = vec3(positions[faces[i][index]]);
            positions.push_back(v);

            // add the corrected tex-coord
            texCoords.push_back(vec2(texExtra));

            // and of course normal needs to be copied
            v = vec3(normals[faces[i][index]]);
            normals.push_back(v);

            // we change the vertex of the face to the newly added vertex, which is at end of positions
            faces[i][index] = positions.size() - 1;

        }

        // check v
        small1 = uv1[1] < 0.25;
        small2 = uv2[1] < 0.25;
        small3 = uv3[1] < 0.25;

        big1 = uv1[1] > 0.75;
        big2 = uv2[1] > 0.75;
        big3 = uv3[1] > 0.75;

        index = -1;
        if(small1 && small2 && big3){
            index = 2;
            texExtra = vec2(uv3[0], uv3[1]-1.0);
        }else if(small1 && big2 && small3){
            index = 1;
            texExtra = vec2(uv2[0], uv2[1]-1.0);
        }else if(big1 && small2 && small3){
            index = 0;
            texExtra = vec2(uv1[0], uv1[1]-1.0);
        }else if(big1 && big2 && small3){
            index = 2;
            texExtra = vec2(uv3[0], uv3[1]+1.0);
        }else if(big1 && small2 && big3){
            index = 1;
            texExtra = vec2(uv2[0], uv2[1]+1.0);
        }else if(small1 && big2 && big3){
            index = 0;
            texExtra = vec2(uv1[0], uv1[1]+1.0);
        }

        if(index != -1){
            positions.push_back(vec3(positions[faces[i][index]]));
            normals.push_back(vec3(normals[faces[i][index]]));
            texCoords.push_back(vec2(texExtra));
            faces[i][index] = positions.size()-1;
            buffered = false;
        }
    }
}

vector<vec3> TriMesh::createHighlightBox(void){
    vector<vec3> highlightBox = vector<vec3>();
    // first face
    // first line
    highlightBox.push_back(vec3(boundingBoxMin));
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMin.y, boundingBoxMax.z));
    // 2nd line
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMin.y, boundingBoxMax.z));
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMax.y, boundingBoxMax.z));
    // 3rd line
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMax.y, boundingBoxMax.z));
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMax.y, boundingBoxMin.z));
    // 4th line
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMax.y, boundingBoxMin.z));
    highlightBox.push_back(vec3(boundingBoxMin));

    // second face
    // first line is equal to first line from last
    // 2nd line
    highlightBox.push_back(vec3(boundingBoxMin));
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMin.y, boundingBoxMin.z));
    // 3rd line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMin.y, boundingBoxMin.z));
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMax.y, boundingBoxMin.z));
    // 4th line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMax.y, boundingBoxMin.z));
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMax.y, boundingBoxMin.z));

    // third face
    // 1st line already done
    // 2nd line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMin.y, boundingBoxMin.z));
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMin.y, boundingBoxMax.z));
    // 3rd line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMin.y, boundingBoxMax.z));
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMax.y, boundingBoxMax.z));
    // 4th line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMax.y, boundingBoxMax.z));
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMax.y, boundingBoxMin.z));

    // we have now described 3 faces of a cube. We only need to create two more lines to connect the open ends of the cube. Everything else is done
    // first line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMax.y, boundingBoxMax.z));
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMax.y, boundingBoxMax.z));
    // second line
    highlightBox.push_back(vec3(boundingBoxMax.x, boundingBoxMin.y, boundingBoxMax.z));
    highlightBox.push_back(vec3(boundingBoxMin.x, boundingBoxMin.y, boundingBoxMax.z));

    return highlightBox;
}

bool TriMesh::intersectTriangle(Ray & ray, mat4 & modelMatrix, vec3 & cameraPosition, vec3 & intersectionPoint, vec3 & barycentricPoint, int & faceIndex){


    // get data from ray
    vec3 origin;
    vec3 direction;
    ray.getOrigin(origin);
    ray.getDirection(direction);

    // first we check for a bounding sphere around our object. if there is no hit, we did not hit any primitive within the mesh
    vec3 center = vec3(0);
    vec3 upVector = vec3(0,0,1);
    center = vec3(modelMatrix * vec4(center,1));
    upVector = vec3(modelMatrix * vec4(upVector,1));
    vec3 inters;
    vec3 norm;

    if(!glm::intersectRaySphere(origin, glm::normalize(direction), center, glm::distance(center, upVector), inters, norm)){
        return false;
    }


    // now transform the vertices into world-coordinates (as light and camera position are in world coordinates)
    vector<vec3> worldPositions = vector<vec3>(positions.size());
    for(unsigned int i=0; i<positions.size(); i++){
        worldPositions[i] = vec3(modelMatrix*vec4(positions[i],1));
    }

    // buffer for found points
    vec3 intersection;
    vec3 barycentricIntersection;


    // found is return value. tells if any triangle found
    bool found = false;
    // foundThis is specific to current triangle
    bool foundThis = false;
    // currently found minimum distance
    float minDist = 10000000.0;
    for(unsigned int i=0; i<faces.size(); i++){
        vec3 v1 = worldPositions[faces[i][0]];
        vec3 v2 = worldPositions[faces[i][1]];
        vec3 v3 = worldPositions[faces[i][2]];

        // with this do-while-loop we ensure that we don't hit ourself on reflection
        bool samePoint = false;
        do{
            foundThis = glm::intersectRayTriangle(origin, direction, v1, v2, v3, barycentricIntersection);

            if(foundThis){
                // intersection position as adapted from https://github.com/g-truc/glm/issues/6
                intersection = origin + direction * barycentricIntersection.z;

                samePoint = glm::distance(origin, intersection) < 0.01;

                origin = origin + 0.1f*direction;
            }

        }while(foundThis && samePoint);

        if(foundThis){
           float distance = glm::abs(glm::distance(cameraPosition, intersection));

           if(found){
               // only replace if this intersection is closer to camera than previous             
               if(distance < minDist){
                   intersectionPoint = vec3(intersection);
                   barycentricPoint = vec3(1.0f-barycentricIntersection.x-barycentricIntersection.y, barycentricIntersection.x, barycentricIntersection.y);
                   minDist = distance;
                   faceIndex = i;
               }
           }else{
               intersectionPoint = vec3(intersection);
               barycentricPoint = vec3(1.0f-barycentricIntersection.x-barycentricIntersection.y, barycentricIntersection.x, barycentricIntersection.y);
               minDist = distance;
               faceIndex = i;
           }
        }

        found = found | foundThis;
    }

    return found;
}

vec3 TriMesh::computeTriangleNormal(vec3 & baryCoords, int &faceIndex){

    // Normalvectors of vertices v1,v2,v3
    vec3 nv1 = normals[faces[faceIndex][0]];
    vec3 nv2 = normals[faces[faceIndex][1]];
    vec3 nv3 = normals[faces[faceIndex][2]];

    //returns normalized interpolated normal of the triangle
    return vec3(glm::normalize(baryCoords[0]*nv1 + baryCoords[1]*nv2 + baryCoords[2]*nv3));
}

void TriMesh::getUVCoords(vec3 & baryCoords, int faceIndex, vec2 & uv){

    vec2 u0 = texCoords[faces[faceIndex][0]];
    vec2 u1 = texCoords[faces[faceIndex][1]];
    vec2 u2 = texCoords[faces[faceIndex][2]];

    //baryCoords = glm::normalize(baryCoords);
    uv = vec2(baryCoords[0] * u0 + baryCoords[1] * u1 + baryCoords[2] * u2);
}

void TriMesh::deleteVAO(){


    if (vaos[0] != 0){
        glDeleteVertexArrays(2,&vaos[0]);
    }

    buffered = false;
}

// generate buffers for given data
void TriMesh::generateVBOBuffers(void){

    vbos = vector<GLuint>(5);

    glGenVertexArrays(2, &vaos[0]);
    glBindVertexArray(vaos[0]);


    /*****************************
    *  Normal Model
    ******************************/


    // vertex buffer
    glGenBuffers(1, &vbos[0]);
    glBindBuffer(GL_ARRAY_BUFFER, vbos[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*this->positions.size(), &this->positions[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(attribVertex);
    glVertexAttribPointer(attribVertex, 3, GL_FLOAT, GL_FALSE, 0, 0);


    glGenBuffers(1, &vbos[1]);
    glBindBuffer(GL_ARRAY_BUFFER, vbos[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*this->normals.size(), &this->normals[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(attribNormal);
    glVertexAttribPointer(attribNormal, 3, GL_FLOAT, GL_FALSE, 0, 0);

    //upload texture coordinates
    //TODO make this conditional
    glGenBuffers(1, &vbos[2]);
    glBindBuffer(GL_ARRAY_BUFFER, vbos[2]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2)*this->texCoords.size(), &this->texCoords[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(attribTexCoord);
    glVertexAttribPointer(attribTexCoord, 2, GL_FLOAT, GL_FALSE, 0, 0);



    // index buffer
    glGenBuffers(1, &vbos[3]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos[3]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(glm::vec3)*this->faces.size(), &this->faces[0], GL_STATIC_DRAW);


    //glBindVertexArray(0);


    /***********************
     *      Highlighting
     ***********************/

    vector<vec3> highlightBox = createHighlightBox();
    glBindVertexArray(vaos[1]);

    glGenBuffers(1, &vbos[4]);
    glBindBuffer(GL_ARRAY_BUFFER, vbos[4]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*highlightBox.size(), &highlightBox[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(attribVertex);
    glVertexAttribPointer(attribVertex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindVertexArray(0);
}



// draw the mesh using vertex arrays
void TriMesh::draw(void){

    if(!buffered){
        generateVBOBuffers();
        buffered = true;
    }

    glBindVertexArray(vaos[0]);
    glDrawElements(GL_TRIANGLES, 3 * faces.size(), GL_UNSIGNED_INT, 0);

    if(selected) {
        glBindVertexArray(vaos[1]);
        glDrawArrays(GL_LINES, 0, 12 * 3);
    }

    glBindVertexArray(0);
}

void TriMesh::generateDefaultUV(){

    texCoords.clear();

    texCoords.push_back(vec2(0, 0));
    texCoords.push_back(vec2(0, 1));
    texCoords.push_back(vec2(1, 1));
    texCoords.push_back(vec2(1, 0));

}

