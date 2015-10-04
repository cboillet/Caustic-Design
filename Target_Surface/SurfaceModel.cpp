#include "global.h"
#include <math.h>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/* io */
#include <iostream>
#include <fstream>


/*local*/
#include "SurfaceModel.h"
#include "SurfaceMesh.h"


Model::Model(GLchar* path){
    //this->loadModel(path);
    this->focalLength = 40;
    this->surfaceSize = 5;
}

Model::Model()
{
    this->focalLength = 40;
    this->surfaceSize = 5;
}

void Model::exportModel(std::string filename)
{

    std::cout << "current number of meshes " << scene->mNumMeshes << std::endl;
    //scene->mNumMeshes = 0;
    //scene->mMeshes[0] = NULL;
    scene->mMeshes = new aiMesh*[1];
    scene->mMeshes[0] = new aiMesh();

    aiMesh* mesh = scene->mMeshes[0];

    // --- prepare scene for exportation. (re-write vertices)
    // prepare mesh
    std::vector<glm::uvec3> faces = meshes[0].indices;
    uint nFaces = faces.size();
    uint nIndices = 3;
    //aiMesh* mesh = new aiMesh();

    std::vector<Vertex> verts = meshes[0].vertices;
    //meshes[0].expandVertices(verts);

    int maxAxis;
    float max = 0.0f;
    for (uint i=0; i<verts.size(); i++)
    {
        for (uint j=0; j<3; j++)
        {
            if(fabs(verts[i].Position[j]) > max)
            {
                max = fabs(verts[i].Position[j]);
                maxAxis = j;
            }
        }
    }

    std::cout << "max = " << max << ", on axis " << maxAxis << std::endl;

    uint nVertices = verts.size();
    mesh->mVertices = new aiVector3D[nVertices];
    mesh->mNormals = new aiVector3D[nVertices];
    mesh->mNumVertices = nVertices;

    for (uint i=0; i<nVertices; i++)
    {
        glm::vec3 pos = verts[i].Position;// * scaling;
        mesh->mVertices[i] = aiVector3D(pos.x/surfaceSize, pos.y/surfaceSize, pos.z/surfaceSize);//verts[i].Position.x, verts[i].Position.y, verts[i].Position.z);
        mesh->mNormals[i] = aiVector3D(verts[i].Normal.x, verts[i].Normal.y, verts[i].Normal.z);
    }

    //mesh->mVertices = vertices;
    //mesh->mNormals = normals;



    mesh->mFaces = new aiFace[ nFaces ];
    mesh->mNumFaces = nFaces;

    for (uint i=0; i<nFaces; i++)
    {
        aiFace& face = mesh->mFaces[i];
        face.mIndices = new uint[nIndices];
        face.mNumIndices = nIndices;

        for (uint j=0; j<nIndices; j++)
        {
            face.mIndices[j] = faces[i][j];
        }
    }

    // --- actually export
    Assimp::Exporter exporter;
    exporter.Export(scene, "stl", filename);

    std::cout << "exported to " << filename << std::endl;

    //for(uint i=0; i<verts.size(); i++)
    //    delete verts[i];

    //delete[] aiMeshes;
    //delete[] vertices;
    //delete[] normals;
    //delete mesh;
}

void Model::loadReceiverLightPoints(QString path)
{
    if(meshes.empty())
        std::cerr << "Load mesh first" << std::endl;

    receiverLightPositions.clear();

    // load 2D Point as 3D point (focal_length is x-coordinate)
    std::ifstream ifs(qPrintable(path));
    float y,z;

    while(ifs >> y >> z) receiverLightPositions.push_back(glm::vec3(focalLength + meshes[0].getMaxX(), y*surfaceSize/CAUSTIC_DOMAIN, z*surfaceSize/CAUSTIC_DOMAIN));

    std::cout << "Loaded " << receiverLightPositions.size() << " light positions for focal length " << focalLength << std::endl;

    //calculate and load the desired normals
    computeLightDirectionsScreenSurface();
    fresnelMapping();
}

void Model::loadModel(string path){
    clean();

    Assimp::Importer import;
    const aiScene* scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_JoinIdenticalVertices);
    if(!scene || scene->mFlags == AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    {
        std::cout << "ERROR::ASSIMP::" << import.GetErrorString() << std::endl;
        return;
    }
    this->directory = path.substr(0, path.find_last_of('/'));

    this->processNode(scene->mRootNode, scene);

    // copy the scene to a modifiable copy.
    aiCopyScene(scene, &this->scene);
    meshes[0].setVertices();
    meshes[0].createEdgeIndices();
    meshes[0].calcEdgeAdjacentFaces();
    meshes[0].createFrontFaceMatrix();
}

void Model::processNode(aiNode *node, const aiScene *scene){
    // Process all the node's meshes (if any)
    for(GLuint i = 0; i < node->mNumMeshes; i++)
    {
        aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
        this->meshes.push_back(this->processMesh(mesh, scene));
    }
    // Then do the same for each of its children
    for(GLuint i = 0; i < node->mNumChildren; i++)
    {
        this->processNode(node->mChildren[i], scene);
    }
}

Mesh Model::processMesh(aiMesh* mesh, const aiScene* scene){
    vector<Vertex> vertices;
    vector<glm::uvec3> indices;
    vector<Texture> textures;
    vector<int> edgeIndicesToIndices;
    float max_X = -100;
    float max_Y = -100;
    float max_Z = -100;
    int edgeIndices;
    int j;
    std::cout << "loading mesh with " << mesh->mNumVertices << " vertices" << std::endl;

    for(GLuint i = 0; i < mesh->mNumVertices; i++)
    {
        Vertex vertex;
        // Process vertex positions, normals and texture coordinates
        glm::vec3 vector;
        vector.x = mesh->mVertices[i].x;
        vector.y = mesh->mVertices[i].y;
        vector.z = mesh->mVertices[i].z;
        vector *= surfaceSize;
        vertex.Position = vector;
        if (max_X<vector.x) max_X=vector.x;
        if (max_Y<vector.y) max_Y=vector.y;
        if (max_Z<vector.z) max_Z=vector.z;
        //normals
        //vector.x = mesh->mNormals[i].x;
        //vector.y = mesh->mNormals[i].y;
        //vector.z = mesh->mNormals[i].z;
        //vertex.Normal = vector;

        if(mesh->mTextureCoords[0]) // Does the mesh contain texture coordinates?
        {
            glm::vec2 vec;
            vec.x = mesh->mTextureCoords[0][i].x;
            vec.y = mesh->mTextureCoords[0][i].y;
            vertex.TexCoords = vec;
        }
        else
            vertex.TexCoords = glm::vec2(0.0f, 0.0f);

        vertices.push_back(vertex);
    }
    // Process indices
    for(GLuint i = 0; i < mesh->mNumFaces; i++)
    {
        aiFace face = mesh->mFaces[i];
        if(face.mNumIndices != 3)
        {
            std::cerr << "can't import model, indices wrong" << std::endl;
        }

        indices.push_back(glm::uvec3(face.mIndices[0], face.mIndices[1], face.mIndices[2]));        
//        if(isEdge(mesh->mVertices[face.mIndices[0]].x) || isEdge(mesh->mVertices[face.mIndices[1]].y) || isEdge(mesh->mVertices[face.mIndices[0]].z) && mesh->mVertices[face.mIndices[1]].x==max_X && mesh->mVertices[face.mIndices[1]].x==max_X && mesh->mVertices[face.mIndices[2]].x==max_X){ //one of the vertices
//            edgeIndicesToIndices.push_back(j);
//            edgeIndices++;
//        }
        //for(GLuint j = 0; j < face.mNumIndices; j++)
        //    indices.push_back(face.mIndices[j]);
        j++;
    }

    // Process material
    if(mesh->mMaterialIndex >= 0)
    {
        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
        vector<Texture> diffuseMaps = this->loadMaterialTextures(material,
                                            aiTextureType_DIFFUSE, "texture_diffuse");
        textures.insert(textures.end(), diffuseMaps.begin(), diffuseMaps.end());
        vector<Texture> specularMaps = this->loadMaterialTextures(material,
                                            aiTextureType_SPECULAR, "texture_specular");
        textures.insert(textures.end(), specularMaps.begin(), specularMaps.end());
    }
    return Mesh(vertices, textures, indices, max_X, max_Y, max_Z/*, edgeIndicesToIndices*/);
}

void Model::rescaleMeshes(float newScale)
{

    // rescale mesh itself
    for (uint i=0; i<meshes.size(); i++)
    {
        for (uint j=0; j<meshes[i].vertices.size(); j++)
        {
            /*glm::vec3 buf = meshes[i].vertices[j].Position;
            buf /= oldScale;
            buf *= newScale;
            meshes[i].vertices[j].Position = buf;*/
            meshes[i].vertices[j].Position = (meshes[i].vertices[j].Position / surfaceSize) * newScale;
        }
        meshes[i].calcMax();
    }

    // rescale light-ray-pos

    for (uint i=0; i<receiverLightPositions.size(); i++)
    {
        receiverLightPositions[i] = glm::vec3(receiverLightPositions[i].x,
                                            newScale*receiverLightPositions[i].y / surfaceSize,
                                            newScale*receiverLightPositions[i].z / surfaceSize);
    }
}

void Model::modifyMesh()
{
    if(meshes.empty()) return;

    int modifyIndex = meshes[0].vertices.size()/3;
    glm::vec3 oldVal =  meshes[0].vertices[modifyIndex].Position;
    meshes[0].vertices[modifyIndex].Position = 2.0f * meshes[0].vertices[modifyIndex].Position;
    glm::vec3 newVal = meshes[0].vertices[modifyIndex].Position;

    std::cout << "modified vertex " << modifyIndex << " from (" <<
                 oldVal.x << ", " << oldVal.y << ", " << oldVal.z  << ") to ("  <<
              newVal.x << ", " << newVal.y << ", " << newVal.z  << ")" << std::endl;

    // update normals
    meshes[0].calculateVertexNormals();
}

void Model::shootRay(vector<glm::highp_dvec3> & direction, vector<glm::highp_dvec3> & redirect, vector<glm::highp_dvec3> & endpoint){
    double refraction = MATERIAL_REFRACTIV_INDEX;
    std::cout<<"faceVerticesEdge size"<<meshes[0].faceVerticesEdge.size()<<std::endl;
    //uint i = 0;
    for(uint i=0; i<meshes[0].faceVerticesEdge.size(); i++){
        //i=in*50;
        glm::highp_dvec3 d = glm::highp_dvec3(1,0,0);
        direction.push_back(d);
        glm::highp_dvec3 redi =  glm::highp_dvec3(meshes[0].vertices[i].Position);
        redirect.push_back(redi);
        glm::highp_dvec3 pointNormal = glm::highp_dvec3(meshes[0].vertices[i].Normal);
        double c = glm::dot(pointNormal, d);
        glm::highp_dvec3 refractedRay = refraction * d + (refraction*c - sqrt(1 - refraction*refraction*(1-c*c))) * (-pointNormal);
        refractedRay = glm::normalize(refractedRay);

        //std::cerr << refractedRay.x << ", " << refractedRay.y << ", " << refractedRay.z << std::endl;

        glm::highp_dvec3 ed = double(1.5*getFocalLength()) * refractedRay + redi;
        endpoint.push_back(ed);
        //std::cerr << endpoint.x << ", " << endpoint.y << ", " << endpoint.z << std::endl;
    }
}

void Model::setFocalLength(float newLength)
{
    this->focalLength = newLength;

    for (uint i=0; i<receiverLightPositions.size(); i++)
    {
        receiverLightPositions[i].x = newLength + meshes[0].getMaxX();
    }
}

void loadToSurface(int index){
//    double m_dx, m_dy;
//    bool ok = m_image.load(filename);
//    //bool ok = m_image.load(filename.toStdString());
//    if (!ok) return false;

//    this->filename = filename;

//    unsigned w = get_width();
//    unsigned h = get_height();
//    m_dx = 0.5;
//    m_dy = m_dx * double(h) / double(w);
}

vector<Texture> Model::loadMaterialTextures(aiMaterial *mat, aiTextureType type, string typeName){
    vector<Texture> textures;
    for(GLuint i = 0; i < mat->GetTextureCount(type); i++)
    {
        aiString str;
        mat->GetTexture(type, i, &str);
        Texture texture;
        texture.id = TextureFromFile(str.C_Str(), this->directory);
        texture.type = typeName;
        //texture.path = str;
        textures.push_back(texture);
    }
    return textures;
}

GLint TextureFromFile(const char* path, string directory)
{
//     //Generate texture ID and load texture data
//    string filename = string(path);
//    filename = directory + '/' + filename;
//    GLuint textureID;
//    glGenTextures(1, &textureID);
//    int width,height;
//    unsigned char* image = SOIL_load_image(filename.c_str(), &width, &height, 0, SOIL_LOAD_RGB);
//    // Assign texture to ID
//    glBindTexture(GL_TEXTURE_2D, textureID);
//    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
//    glGenerateMipmap(GL_TEXTURE_2D);

//    // Parameters
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//    glBindTexture(GL_TEXTURE_2D, 0);
//    SOIL_free_image_data(image);
//    return textureID;
}


void Model::printAllVertices(){
    for (int j=0; j<36; j++){
        if(j%3 == 0) std::cout<<"\n"<<std::endl;
        std::cout<<"vertice"<<j<<" x:"<<meshes[0].vertices[j].Position.x<<" y:"<<meshes[0].vertices[j].Position.y<<" z:"<<meshes[0].vertices[j].Position.z<<std::endl;
        }
    //meshes[0].getSpatialLayout(0);
}

int Model::findSurfaceMesh(){
    int area;
    double temp;
    for (int i=0; i<meshes.size(); i++){
  //      meshes[i].getSurface();
    }
}


void Model::setNormals(bool edge){
    int limit;
    currentNormals.clear();
    if (edge) limit = meshes[0].faceVerticesEdge.size();
    else limit = meshes[0].faceVerticesEdge.size();
    for (int i=0; i<meshes[0].faceVerticesEdge.size(); i++){
        glm::vec3 norm ;
        norm.x = meshes[0].faceVerticesEdge[i]->Normal.x;
        norm.y = meshes[0].faceVerticesEdge[i]->Normal.y;
        norm.z = meshes[0].faceVerticesEdge[i]->Normal.z;

        currentNormals.push_back(norm);
    }
}

void Model::computeLightDirectionsScreenSurface(){
    //get the position on the surface without the corner vertices
    int j=0;
    screenDirections.clear();
    glm::vec3 vecNorm;
    for(int i=0; i<meshes[0].faceVertices.size(); i++){
        //compute vector surface to screen
//         if(!meshes[0].isEdge(meshes[0].faceVerticesEdge[i])){
//            vecNorm = receiverLightPositions[j] - meshes[0].faceVerticesEdge[i]->Position ;
//            j++;
//            std::cout<<"load point: "<<i<<"from light receiver position"<<j<<std::endl;
//         }
//         else {
//             vecNorm = meshes[0].faceVerticesEdge[i]->Normal;
//             std::cout<<"load edge: "<<i<<std::endl;
//         }
        vecNorm = receiverLightPositions[i] - meshes[0].faceVertices[i]->Position ;
        screenDirections.push_back(glm::normalize(vecNorm));
    }
}

//compute the desired normals
void Model::fresnelMapping(){
    //calculate sin(i1)/sin(i2)
    float refraction = MATERIAL_REFRACTIV_INDEX;
    desiredNormals.clear();
    int j;
    for(int i = 0; i<meshes[0].faceVertices.size(); i++){
        glm::vec3 incidentLight;
        incidentLight.x = -1;
        incidentLight.y = 0;
        incidentLight.z = 0;

        glm::vec3 vert = screenDirections[i]/refraction;
        //vert *= refraction;
        glm::vec3 norm;

        //normal of the surface see Kiser and Pauly
        //norm = (incidentLight +  (1/refraction)*vert)/(glm::length(incidentLight+(1/refraction)*vert));
        norm = (incidentLight + vert);
        //norm = (incidentLight + refraction*vert);
        desiredNormals.push_back(glm::normalize(norm)*-1.0f);
        //desiredNormals.push_back(glm::normalize(norm)*-1.0f);


    }
}
