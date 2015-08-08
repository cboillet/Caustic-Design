#include "global.h"
#include <math.h>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*local*/
#include "SurfaceMesh.h"
#include "utils.h"



Mesh::Mesh(vector<Vertex> vertices, vector<Texture> textures)
{
    this->vertices = vertices;
    this->textures = textures;
    create_indices();
    shrink_vertices();
    calcMaxX();
}

void Mesh::create_indices()
{
    // first we create an array with a 1:1 mapping
    uint faces = vertices.size() / 3;

    indices.resize(faces);
    for (uint i=0; i<faces; i++)
    {
        uint base = 3*i;
        indices[i] = glm::uvec3(base, base+1, base+2);
    }

}


void Mesh::shrink_vertices()
{
    // we (for each vertex in each face) check if there is a vertex with a lower index that we can take instead. And we keep track of maximum index used and what values are used
    bool* used=new bool[vertices.size()];

    for (uint i=0; i<vertices.size(); i++)
    {
        used[i] = false;
    }

    for (int i=vertices.size()-1; i >= 0; i--)
    {
        for (int j=0; j<i; j++)
        {

            // if distance between the two vertices is small, we assume they are the same -> set higher index to lower one
            if(glm::distance(vertices[i].Position, vertices[j].Position) < 0.00001){
                int indexBase = i/3;
                int indexOffset = i%3;

                indices[indexBase][indexOffset] = j;
                used[j] = true;

                break;
            }
        }
    }

    // we now always point to the vertex with the lowest index
    // and we know the highest index we point to in the vertex-vector.
    // So we can remove everything that is unused

    // WARNING: Each index that is higher than the removed index needs to be decreased by one
    // we start at max index, which is important since we then don't need to update the used[] array
     for (int i=vertices.size()-1; i>=0; i--)
    {
        if(!used[i])
        {
            // remove from vertex-vector and decrease indices
            vertices.erase(vertices.begin()+i);

            // decrease if needed
            for (uint index=0; index<indices.size(); index++)
            {
                // we use uvec3, so lenght is always 3
                for (uint vIndex=0; vIndex<3; vIndex++)
                {
                    // no need to check for equality, since we know the vertex at current index is not used
                    if(indices[index][vIndex] > i)
                        indices[index][vIndex] = indices[index][vIndex]-1;
                }
            }
        }
    }

    std::cout << "vertices.size reduced to " << vertices.size() << std::endl;

    delete[] used;
}

void Mesh::calcMaxX()
{
    maxX = -std::numeric_limits<float>::max();

    foreach (Vertex v, vertices)
    {
        if(v.Position.x > maxX)
            maxX = v.Position.x;
    }

    std::cout << "maxX = " << maxX << std::endl;
}


float Mesh::vectorNorm(Vertex v1, Vertex v2){
    int sum = pow((v2.Position.z-v1.Position.z),2)+pow((v2.Position.y-v1.Position.y),2)+pow((v2.Position.x-v1.Position.x),2);
    return sqrt(sum);
}

void Mesh::rescale(float oldScale, float newScale)
{
    for (uint i=0; i<vertices.size(); i++)
    {
        //vertices[i].Position = (vertices[i].Position)
    }
}

vector<int> Mesh::longerSegment(Vertex v0, Vertex v1, Vertex v2, int first){
    vector<int> ind;
    int max = vectorNorm(v0,v1);
    if(vectorNorm(v1,v2)>max) {
        ind.push_back(first + 1);
        ind.push_back(first + 2);
        ind.push_back(first); //no concerned by longest segment
    }
    else if(vectorNorm(v0,v2)>max) {
        ind.push_back(first);
        ind.push_back(first + 2);
        ind.push_back(first + 1);
    }
    else {
        ind.push_back(first);
        ind.push_back(first + 1);
        ind.push_back(first + 2);
    }
    return ind;
}

bool Mesh::adjacent(vector<Vertex> vec1, vector<Vertex> vec2){
    bool adj = false;
    int equal;
    for(int i=0; i<3; i++){
        for (int j=0; j<3; j++)
            if ((vec1[i].Position.x == vec2[j].Position.x) && (vec1[i].Position.y == vec2[j].Position.y) && (vec1[i].Position.z == vec2[j].Position.z)) equal++;
    }
    if (equal >= 2) adj = true;
    return adj;
}

bool Mesh::equal(vector<Vertex> vec1, vector<Vertex> vec2){
    bool adj = false;
    int equal;
    for(int i=0; i<3; i++){
        for (int j=0; j<3; j++)
            if ((vec1[i].Position.x == vec2[j].Position.x) && (vec1[i].Position.y == vec2[j].Position.y) && (vec1[i].Position.z == vec2[j].Position.z)) equal++;
    }
    if (equal >= 3) adj = true;
    return adj;
}

bool Mesh::compareArea(vector<Vertex> vec1, vector<Vertex> vec2){
    int a1,a2;
    vector<int> ia1 = longerSegment(vec1[0],vec1[1],vec1[2],0);
    vector<int> ia2 = longerSegment(vec2[0],vec2[1],vec2[2],0);
    int edgea1 = outTriplet(ia1,0,2);
    int edgea2 = outTriplet(ia2,0,2);
    a1 = vectorNorm(vec1[ia1[0]],vec1[edgea1])*vectorNorm(vec1[ia1[1]],vec1[edgea1])/2;
    a2 = vectorNorm(vec2[ia2[0]],vec2[edgea2])*vectorNorm(vec2[ia2[1]],vec2[edgea2])/2;
    if (a2>=a1) return true;
    else return false;
}

void Mesh::exportVertices(const QString& filename){
    std::ofstream ofs(qPrintable(filename));
    ofs.precision(20);
    vector<Vertex> verticesMesh = selectVerticesMeshFace();
    int exported = 0, excluded = 0;
    for (unsigned i = 0; i < verticesMesh.size(); ++i)
    {
        // don't export vertices on the edge
        ofs << verticesMesh[i].Position.y << " " << verticesMesh[i].Position.z << std::endl;

    }
    ofs.close();
}


vector<Vertex>  Mesh::selectVerticesMeshFace(){
    vector<Vertex> faceVertex;
    int max =  0;
    for (int i = 0; i < vertices.size(); ++i){
        if (vertices[i].Position.x>max) max = vertices[i].Position.x;
    }

    for (int i = 0; i < vertices.size(); ++i){
         if(abs(vertices[i].Position.x-max) < 0.00001){
             // we are on side with max x-values. now exclude min/max y- and z- values
             if(!floatEquals(fabs(vertices[i].Position.y), 1.0f) && !floatEquals(fabs(vertices[i].Position.z), 1.0f))
             faceVertex.push_back(vertices [i]);
         }
    }
    return faceVertex;
}

