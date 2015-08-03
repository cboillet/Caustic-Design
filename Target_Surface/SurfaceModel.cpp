#include "global.h"
#include <math.h>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*openGL*/
#include "SOIL.h"

/*local*/
#include "SurfaceModel.h"
#include "SurfaceMesh.h"


Model::Model(GLchar* path){
    this->loadModel(path);
}

void Model::loadModel(string path){
    Assimp::Importer import;
    const aiScene* scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);

    if(!scene || scene->mFlags == AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    {
        std::cout << "ERROR::ASSIMP::" << import.GetErrorString() << std::endl;
        return;
    }
    this->directory = path.substr(0, path.find_last_of('/'));

    this->processNode(scene->mRootNode, scene);

    //this->Draw();
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
    vector<GLuint> indices;
    vector<Texture> textures;

    for(GLuint i = 0; i < mesh->mNumVertices; i++)
    {
        Vertex vertex;
        // Process vertex positions, normals and texture coordinates
        glm::vec3 vector;
        vector.x = mesh->mVertices[i].x;
        vector.y = mesh->mVertices[i].y;
        vector.z = mesh->mVertices[i].z;
        vertex.Position = vector;
        //normals
        vector.x = mesh->mNormals[i].x;
        vector.y = mesh->mNormals[i].y;
        vector.z = mesh->mNormals[i].z;
        vertex.Normal = vector;

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
        for(GLuint j = 0; j < face.mNumIndices; j++)
            indices.push_back(face.mIndices[j]);
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
    return Mesh(vertices, indices, textures);
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
    meshes[0].getSpatialLayout(0);
}

int Model::findSurfaceMesh(){
    int area;
    double temp;
    for (int i=0; i<meshes.size(); i++){
  //      meshes[i].getSurface();
    }
}



//SurfaceModel::SurfaceModel(int w, int h, int d, string p):width(w),height(h),depth(d),path(p)
//{
        
//}
