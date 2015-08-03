#include "global.h"
#include <math.h>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*local*/
#include "SurfaceMesh.h"

Mesh::Mesh(vector<Vertex> vertices, vector<GLuint> indices, vector<Texture> textures)
{
    this->vertices = vertices;
    this->indices = indices;
    this->textures = textures;
    nbMeshLayout = 2; //nb of triangle mesh in the space layout

    GLenum code = glewInit();
    if(code != GLEW_OK)
    {
        fprintf(stderr, "impossible d'initialiser GLEW : %s\n",
                        glewGetErrorString(code));
    }
}


int Mesh::getSpatialLayout(int nbFace){
//    double ans;
//    Vertex v0 = this->vertices[0];
//    std::cout<<"v0 position x"<<v0.Position.x<<"v0 position y"<<v0.Position.y<<"v0 position z"<<v0.Position.z<<std::endl;
//    int last= vertices.size()-1;
//    Vertex vL = this->vertices[last];
//    std::cout<<' '<<" x:"<<vL.Position.x-v0.Position.x<<std::endl;
//    std::cout<<' '<<" y:"<<vL.Position.y-v0.Position.y<<std::endl;
//    std::cout<<' '<<" z:"<<vL.Position.z-v0.Position.z<<std::endl;
      std::cout<<"Indexes from "<<nbFace*6<<" to "<<(nbFace*6+nbMeshLayout*3)<<std::endl;
      return nbFace*6+nbMeshLayout*3;
}

float Mesh::vectorNorm(Vertex v1, Vertex v2){
    int sum = pow((v2.Position.z-v1.Position.z),2)+pow((v2.Position.y-v1.Position.y),2)+pow((v2.Position.x-v1.Position.x),2);
    return sqrt(sum);
}

vector<int> Mesh::longerSegment(Vertex v0, Vertex v1, Vertex v2){
    vector<int> ind;
    int max = vectorNorm(v0,v1);
    if(vectorNorm(v1,v2)>max) {
        ind.push_back(1);
        ind.push_back(2);
    }
    else if(vectorNorm(v0,v2)>max) {
        ind.push_back(0);
        ind.push_back(2);
    }
    else {
        ind.push_back(0);
        ind.push_back(1);
    }
    return ind;
}

void Mesh::generateTriangles(){
    int edge1, edge2, edge3;
    int index=0;
    int face,test;
    vector<Vertex> vert = vertices;
    std::vector<Vertex>::iterator it;
    it = vert.begin();
    vector<GLuint> ind;
    vector<Texture> text;

    Vertex vertex;
    glm::vec3 vect;

    //1. process vertex position
    vector<int> in = longerSegment(vertices[0],vertices[1],vertices[2]);
    int i1 = 0;
    int i2 = 2;
    std::cout<<"longer segment indice 1"<<i1<<std::endl;
    std::cout<<"longer segment indice 2"<<i2<<std::endl;
    vect.x = (vertices[i1].Position.x+vertices[i2].Position.x)/2;
    vect.y = (vertices[i1].Position.y+vertices[i2].Position.y)/2;
    vect.z = (vertices[i1].Position.z+vertices[i2].Position.z)/2;
    vertex.Position = vect;

    vect.x = (vertices[i1].Normal.x+vertices[i2].Normal.x)/2;
    vect.y = (vertices[i1].Normal.y+vertices[i2].Normal.y)/2;
    vect.z = (vertices[i1].Normal.z+vertices[i2].Normal.z)/2;
    vertex.Normal = vect;

    vertex.TexCoords = glm::vec2(0.0f, 0.0f);
    it = vert.insert(it+(1), vertex);
    it = vert.insert(it+(2), vert[3]);
    it = vert.insert(it+(2), vert[2]);

    //2. process indices
    GLuint lastInd = indices.back();
    indices.push_back((lastInd+1));
//            Texture lastText = textures.back();
//            Texture t2 = lastText;
//            textures.push_back(t2);
    //3. process material
    index+=3;
    test+=1;
    //parseDiffMesh(vert, vertices);
    vertices.clear();
    vertices = vert;

//    while(nbMeshLayout<4)//MESH_AMOUNT && index<vertices.size())
//    {

//        for(int i = 0; i<nbMeshLayout; i++){
//            Vertex vertex;
//            glm::vec3 vect;

//            //1. process vertex position
//            vector<int> in = longerSegment(vertices[index],vertices[index+1],vertices[index+2]);
//            int i1 = 0;
//            int i2 = 2;
//            std::cout<<"longer segment indice 1"<<i1<<std::endl;
//            std::cout<<"longer segment indice 2"<<i2<<std::endl;
//            vect.x = (vertices[i1].Position.x+vertices[i2].Position.x)/2;
//            vect.y = (vertices[i1].Position.y+vertices[i2].Position.y)/2;
//            vect.z = (vertices[i1].Position.z+vertices[i2].Position.z)/2;
//            vertex.Position = vect;

//            vect.x = (vertices[i1].Normal.x+vertices[i2].Normal.x)/2;
//            vect.y = (vertices[i1].Normal.y+vertices[i2].Normal.y)/2;
//            vect.z = (vertices[i1].Normal.z+vertices[i2].Normal.z)/2;
//            vertex.Normal = vect;

//            vertex.TexCoords = glm::vec2(0.0f, 0.0f);
//            it = vert.insert(it+(index+1), vertex);
//            //2. process indices
//            GLuint lastInd = indices.back();
//            indices.push_back((lastInd+1));
////            Texture lastText = textures.back();
////            Texture t2 = lastText;
////            textures.push_back(t2);
//            //3. process material
//            index+=3;
//            test+=1;
//        }
//        vertices.clear();
//        vertices = vert;
//        nbMeshLayout=pow(nbMeshLayout,2);
//        std::cout<<"nb de triangle"<<nbMeshLayout<<std::endl;
//    }

    for (int j=0; j<vertices.size(); j++){
        if(j%3 == 0) std::cout<<"\n"<<std::endl;
        std::cout<<"vertice"<<j<<" x:"<<vertices[j].Position.x<<" y:"<<vertices[j].Position.y<<" z:"<<vertices[j].Position.z<<std::endl;
    }
}

//void Mesh::findEqual(int[] v, vector<Vertex>v2,vector<Vertex>::iterator it){
//    if(it == v2.end()) std::cout<<"x differ"<<int[0]<<int[1]<<int[2]<<std::endl;
//    if(v[0] == (*it).Position.x) && (v[1] == (*it).Position.y) && (v[2] == (*it).Position.z) return;
//    else(findEqual(v), it++);
//}

void Mesh::parseDiffMesh(vector<Vertex> v1, vector<Vertex> v2){
    //a recoder
    int s = max(v1.size(),v2.size());
    int rank=0;
    std::vector<Vertex>::iterator it = v2.begin();
    for (int i=0; i<s; i+=3){
        int[] vertex1 = [v1[i+j].Position.x, v1[i+j].Position.y, v1[i+j].Position.z];
        if (vertex1[0] != v2[rank].Position.x) || (vertex1[1] != v2[rank].Position.y) || (vertex1[0] != v2[rank].Position.z){
     //       findEqual(vertex1, v2,it);
        }
        it++;


//        for(int j=0; j<3; j++){
//            if (v1[i+j].Position.x != v2[i+j].Position.x) {
//                std::cout<<"x differ"<<v1[i+j].Position.x<<v2[i+j].Position.x<<std::endl;
//                bool e = false;
//                equal.push_back(e);
//            }
//            if (v1[i+j].Position.y != v2[i+j].Position.y) {
//                std::cout<<"y differ"<<v1[i+j].Position.y<<v2[i+j].Position.y<<std::endl;
//                bool e = false;
//                equal.push_back(e);
//            }
//            if (v1[i+j].Position.z != v2[i+j].Position.z) {
//                std::cout<<"z differ"<<v1[i+j].Position.z<<v2[i+j].Position.z<<std::endl;
//                bool e = false;
//                equal.push_back(e);
//            }
//            else {
//                bool e = true;
//                equal.push_back(e);
//            }

        }
    }
}

void Mesh::setUpMesh(int nbvertices){
        vector<Vertex> verticesIncluded;
        for (int i=(vertices.size()/6)*5; i<vertices.size();i++){
            Vertex vtemp=vertices[i];
            verticesIncluded.push_back(vtemp);
        }
        //TO DO: improve dimension input to user defined values
        double dx = verticesIncluded[5].Position.x-verticesIncluded[0].Position.x;
        double dy = verticesIncluded[6].Position.x-verticesIncluded[0].Position.x;
        int nbWidth = floor(MESH_AMOUNT) / 3;
        int nbHeight = MESH_AMOUNT / nbWidth;
        double stepx = 2.0 * dx / nbWidth;
        double stepy = 2.0 * dy / nbHeight;

    //    std::vector<Point> points;
    //    for (unsigned i = 0; i < nx; ++i)
    //    {
    //        FT x = (i + 0.5)*stepx - m_domain.get_dx();
    //        x += EPS;
    //        for (unsigned j = 0; j < ny; ++j)
    //        {
    //            FT y = (j + 0.5)*stepy - m_domain.get_dy();
    //            y += EPS;
    //            points.push_back(Point(x, y));
    //        }
    //    }
    //    std::vector<FT> weights(points.size(), 0.0);

    //    for (i = 0; i<nbvertices ; i++)
    //    {
    //        Vertex vertex;
    //        // Process vertex positions, normals and texture coordinates
    //        glm::vec3 vector;
    //        vector.x = mesh->mVertices[i].x;
    //        vector.y = mesh->mVertices[i].y;
    //        vector.z = mesh->mVertices[i].z;
    //        vertex.Position = vector;
    //        //normals
    //        vector.x = mesh->mNormals[i].x;
    //        vector.y = mesh->mNormals[i].y;
    //        vector.z = mesh->mNormals[i].z;
    //        vertex.Normal = vector;

    //        if(mesh->mTextureCoords[0]) // Does the mesh contain texture coordinates?
    //        {
    //            glm::vec2 vec;
    //            vec.x = mesh->mTextureCoords[0][i].x;
    //            vec.y = mesh->mTextureCoords[0][i].y;
    //            vertex.TexCoords = vec;
    //        }
    //        else
    //            vertex.TexCoords = glm::vec2(0.0f, 0.0f);

    //        vertices.push_back(vertex);
    //    }
    //    // Process indices
    //    for(GLuint i = 0; i < mesh->mNumFaces; i++)
    //    {
    //        aiFace face = mesh->mFaces[i];
    //        for(GLuint j = 0; j < face.mNumIndices; j++)
    //            indices.push_back(face.mIndices[j]);
    //    }

    //    // Process material
    //    if(mesh->mMaterialIndex >= 0)
    //    {
    //        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
    //        vector<Texture> diffuseMaps = this->loadMaterialTextures(material,
    //                                            aiTextureType_DIFFUSE, "texture_diffuse");
    //        textures.insert(textures.end(), diffuseMaps.begin(), diffuseMaps.end());
    //        vector<Texture> specularMaps = this->loadMaterialTextures(material,
    //                                            aiTextureType_SPECULAR, "texture_specular");
    //        textures.insert(textures.end(), specularMaps.begin(), specularMaps.end());
    //    }
}
