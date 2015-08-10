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
    calculateVertexNormals();
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

void Mesh::expandVertices(std::vector<Vertex> &outVertices)
{
    for (uint i=0; i<indices.size(); i++)
    {
        // triangles, so index max is 2
        for (uint j=0; j<3; j++)
        {
            outVertices.push_back(vertices[indices[i][j]]);
        }
    }

    std::cout << "expanded vertex-amount from " << vertices.size() << " to " << outVertices.size() << std::endl;
}

void Mesh::calculateFaceNormals(std::vector<glm::vec3> &normals)
{
    for (uint i=0; i<indices.size(); i++)
    {
        glm::vec3 v1, v2;
        v1 = vertices[indices[i][1]].Position - vertices[indices[i][0]].Position;
        v2 = vertices[indices[i][2]].Position - vertices[indices[i][0]].Position;

        normals.push_back(glm::normalize(glm::cross(v1, v2)));
    }
}

void Mesh::calculateVertexNormals()
{
    std::vector<glm::vec3> faceNormals;
    calculateFaceNormals(faceNormals);

    //aF contains indices of adjacend faces per vertex
    vector<vector<unsigned int> > aF;
    aF.resize(vertices.size());
    for(uint i = 0; i < indices.size(); i++){
        aF[indices[i][0]].push_back(i);
        aF[indices[i][1]].push_back(i);
        aF[indices[i][2]].push_back(i);
    }

    //interpolate normals of adjacend faces per vertex
    for(uint i = 0; i< aF.size(); i++){

        glm::vec3 vertexNormal = glm::vec3(0);
        for(uint j= 0; j < aF[i].size(); j++){

            // find out which vertex of the current face is the vertex we are currently looking at
            // aF[i] is a list of faces (aka a list of indices of the indices-vector)
            // so indices[aF[i][j]] is a glm::vec3 that contains one face
            // and the current vertex is vertex[i]
            int thisVertexIndex = -1;
            for(int vIndex=0; vIndex < 3; vIndex++)
            {
                glm::vec3 pos1 = vertices[indices[aF[i][j]][vIndex]].Position;
                glm::vec3 pos2 = vertices[i].Position;
                if(glm::distance(pos1, pos2) < 0.0001f)
                {
                    thisVertexIndex = vIndex;
                    break;
                }

            }

            // we got index of our current vertex within the face, now get others
            int other1 = (thisVertexIndex+1) % 3;
            int other2 = (thisVertexIndex+2) % 3;

            // create the vectors the represent the edges from current vertex to the other 2
            glm::vec3 edge1 = vertices[indices[aF[i][j]][thisVertexIndex]].Position - vertices[indices[aF[i][j]][other1]].Position;
            glm::vec3 edge2 = vertices[indices[aF[i][j]][thisVertexIndex]].Position - vertices[indices[aF[i][j]][other2]].Position;

            // get angle between the edges
            float incidentAngle = abs(glm::angle(glm::normalize(edge1), glm::normalize(edge2)));
            if(incidentAngle > 180)
               incidentAngle = 360 - incidentAngle;

            // use that angle as weighting
            vertexNormal += (faceNormals[aF[i][j]] * incidentAngle);
        }
        vertices[i].Normal = glm::normalize(vertexNormal);
    }
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

void Mesh::exportVertices(const QString& filename, float scaling){
    std::ofstream ofs(qPrintable(filename));
    ofs.precision(20);
    vector<Vertex> verticesMesh = selectVerticesMeshFaceNoEdge();
    int exported = 0, excluded = 0;
    for (unsigned i = 0; i < verticesMesh.size(); ++i)
    {
        // don't export vertices on the edge
        ofs << verticesMesh[i].Position.y * CAUSTIC_DOMAIN / scaling << " " << verticesMesh[i].Position.z  * CAUSTIC_DOMAIN / scaling << std::endl;

    }
    ofs.close();
}


vector<Vertex>  Mesh::selectVerticesMeshFaceNoEdge(){
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

vector<Vertex> Mesh::selectVerticesMeshFaceEdge(){
    vector<Vertex> faceVertex;
    int max =  0;
    for (int i = 0; i < vertices.size(); ++i){
        if (vertices[i].Position.x>max) max = vertices[i].Position.x;
    }

    for (int i = 0; i < vertices.size(); ++i){
         if(abs(vertices[i].Position.x-max) < 0.00001){
             faceVertex.push_back(vertices [i]);
         }
    }
    return faceVertex;
}
