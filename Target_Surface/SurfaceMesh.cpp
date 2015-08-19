#include "global.h"
#include <math.h>
#include <algorithm>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*local*/
#include "SurfaceMesh.h"
#include "utils.h"



Mesh::Mesh(vector<Vertex> vertices, vector<Texture> textures, vector<glm::uvec3> indices)
{
    this->vertices = vertices;
    this->textures = textures;
    this->indices = indices;

    //create_indices();
    //shrink_vertices();
    getAdjacentFacesVector();
    calcMax();
    calculateVertexNormals();
}

void Mesh::create_indices()
{
    std::cout << "creating 1:1 indices.. ";
    std::cout.flush();

    // first we create an array with a 1:1 mapping
    uint faces = vertices.size() / 3;

    indices.resize(faces);
    for (uint i=0; i<faces; i++)
    {
        uint base = 3*i;
        indices[i] = glm::uvec3(base, base+1, base+2);
    }

    //std::cout << "done." << std::endl;
}

std::vector<int> Mesh::createEdgeToNoEdgeMapping()
{
    std::vector<int> mapping;

    for(uint i=0; i<faceVerticesEdge.size(); i++)
    {
        mapping.push_back(getIndexTargetSurface(faceVerticesEdge[i]));
    }

    return mapping;
}

std::vector<int> Mesh::createNoEdgeToEdgeMapping()
{
    std::vector<int> mapping;

    for(uint i=0; i<faceVertices.size(); i++)
    {
        //mapping.push_back(getIndexWithEdge(faceVertices[i]));
    }

    return mapping;
}


void Mesh::shrink_vertices()
{

    std::cout << "shrinking vertices.. ";
    std::cout.flush();

    // we (for each vertex in each face) check if there is a vertex with a lower index that we can take instead. And we keep track of maximum index used and what values are used
    bool* used=new bool[vertices.size()];

    for (uint i=0; i<vertices.size(); i++)
    {
        used[i] = false;
    }

    int percentage = -1;

    for (int i=vertices.size()-1; i >= 0; i--)
    {
        float p = float(vertices.size()-i) / float(vertices.size());
        int currentPercentage = p * 100.0f;
        if(currentPercentage > percentage)
        {
            percentage = currentPercentage;
            std::cout << '\r' << "shrinking vertices.. " << percentage << "% ";
            std::cout.flush();
        }

        //std::cout << '\r' << "shrinking vertices.. " << vertices.size()-i << "/" << vertices.size();

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

    std::cout << "updating indices.. ";
    std::cout.flush();

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
    //std::cout << "done." << std::endl;
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

void Mesh::updateNormal(Vertex *v)
{
    uint vertexIndex;
    std::vector<glm::uvec3> faces;
    std::vector<glm::vec3> normals;
    for (uint i=0; i<indices.size(); i++)
    {

        bool included = false;
        for (uint j=0; j<3; j++)
        {
            if(&vertices[indices[i][j]] == v)
            {
                faces.push_back(indices[i]);
                vertexIndex = indices[i][j];

                glm::vec3 normal;
                calculateFaceNormal(normal, i);
                normals.push_back(normal);

                included = true;
                break;
            }
        }
        if(!included)
        {
            normals.push_back(glm::vec3(0));
        }
    }
    //std::cout << "updating normal " << vertexIndex << std::endl;

    calculateVertexNormal(normals, vertexIndex);
}


void Mesh::calculateFaceNormal(glm::vec3 & normal, uint faceIndex)
{
    glm::vec3 v1, v2;
    v1 = vertices[indices[faceIndex][1]].Position - vertices[indices[faceIndex][0]].Position;
    v2 = vertices[indices[faceIndex][2]].Position - vertices[indices[faceIndex][0]].Position;

    normal = glm::normalize(glm::cross(v1, v2));
}

void Mesh::calculateFaceNormals(std::vector<glm::vec3> &normals)
{
    //std::cout << "calculating face normals.. ";
    //std::cout.flush();

    for (uint i=0; i<indices.size(); i++)
    {
        glm::vec3 normal;
        calculateFaceNormal(normal, i);
        normals.push_back(normal);
    }

   //std::cout << "done." << std::endl;
}

void Mesh::getAdjacentFacesVector()
{
    adjacentFaces.resize(vertices.size());
    for(uint i = 0; i < indices.size(); i++){
        adjacentFaces[indices[i][0]].push_back(i);
        adjacentFaces[indices[i][1]].push_back(i);
        adjacentFaces[indices[i][2]].push_back(i);
    }
}

void Mesh::calcEdgeAdjacentFaces()
{
    edgeAdjacentFaces.resize(edgeIndices.size());

    for(uint i=0; i<edgeIndices.size(); i++)
    {
        edgeAdjacentFaces[edgeIndices[i][0]].push_back(i);
        edgeAdjacentFaces[edgeIndices[i][1]].push_back(i);
        edgeAdjacentFaces[edgeIndices[i][2]].push_back(i);
    }
}

void Mesh::calculateVertexNormal(std::vector<glm::vec3> & faceNormals, uint vertexIndex)
{

    glm::vec3 vertexNormal = glm::vec3(0);
    for(uint j= 0; j < adjacentFaces[vertexIndex].size(); j++){

        // find out which vertex of the current face is the vertex we are currently looking at
        // aF[vertexIndex] is a list of faces (aka a list of indices of the indices-vector)
        // so indices[aF[vertexIndex][j]] is a glm::vec3 that contains one face
        // and the current vertex is vertex[vertexIndex]
        int thisVertexIndex = -1;
        for(int vIndex=0; vIndex < 3; vIndex++)
        {
            glm::vec3 pos1 = vertices[indices[adjacentFaces[vertexIndex][j]][vIndex]].Position;
            glm::vec3 pos2 = vertices[vertexIndex].Position;
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
        glm::vec3 edge1 = vertices[indices[adjacentFaces[vertexIndex][j]][thisVertexIndex]].Position - vertices[indices[adjacentFaces[vertexIndex][j]][other1]].Position;
        glm::vec3 edge2 = vertices[indices[adjacentFaces[vertexIndex][j]][thisVertexIndex]].Position - vertices[indices[adjacentFaces[vertexIndex][j]][other2]].Position;

        // get angle between the edges
        float incidentAngle = abs(glm::angle(glm::normalize(edge1), glm::normalize(edge2)));
        if(incidentAngle > 180)
           incidentAngle = 360 - incidentAngle;

        // use that angle as weighting
        vertexNormal += (faceNormals[adjacentFaces[vertexIndex][j]] * incidentAngle);
    }
    vertices[vertexIndex].Normal = glm::normalize(vertexNormal);
}

void Mesh::calculateVertexNormals()
{
    std::vector<glm::vec3> faceNormals;
    calculateFaceNormals(faceNormals);

    //std::cout << "calculating vertex normals.. ";
    //std::cout.flush();

    //aF contains indices of adjacend faces per vertex
    //interpolate normals of adjacend faces per vertex
    for(uint i = 0; i< vertices.size(); i++){

        calculateVertexNormal(faceNormals, i);
        /*
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
        vertices[i].Normal = glm::normalize(vertexNormal);*/
    }

    //std::cout << "done." << std::endl;
}


void Mesh::calcMax()
{

    std::cout << "calculating max.. ";

    maxX = -std::numeric_limits<float>::max();
    maxY = -std::numeric_limits<float>::max();
    maxZ = -std::numeric_limits<float>::max();

    foreach (Vertex v, vertices)
    {
        if(abs(v.Position.x) > maxX)
            maxX = abs(v.Position.x);

        if(abs(v.Position.y) > maxY)
            maxY = abs(v.Position.y);

        if(abs(v.Position.z) > maxZ)
            maxZ = abs(v.Position.z);
    }

    //std::cout << "done" << std::endl;
    std::cout << "maxX = " << maxX << ", maxX = " << maxY << ", maxZ = " << maxZ << std::endl;
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
    vector<Vertex*> verticesMesh = faceVertices;
    int exported = 0, excluded = 0;
    for (unsigned i = 0; i < verticesMesh.size(); ++i)
    {
        // don't export vertices on the edge
        ofs << verticesMesh[i]->Position.y * CAUSTIC_DOMAIN / scaling << " " << verticesMesh[i]->Position.z  * CAUSTIC_DOMAIN / scaling << std::endl;

    }
    ofs.close();
}

void printVertex(Vertex v)
{
    std::cout << "[" << v.Position.x << ", " << v.Position.y << ", " << v.Position.z << "]";
}

bool Mesh::isEdge(Vertex * v){
    bool edge = false;
    if(floatEquals(fabs(v->Position.y), maxY) || floatEquals(fabs(v->Position.z), maxZ)) edge=true;
    return edge;

}

vector<Vertex*>  Mesh::selectVerticesMeshFaceNoEdge(){

    vector<Vertex*> retVal;

    for(uint i=0; i<vertices.size(); i++)
    {
        if(floatEquals(vertices[i].Position.x, maxX) && !floatEquals(fabs(vertices[i].Position.y), maxY) && !floatEquals(fabs(vertices[i].Position.z), maxZ))
        {
            retVal.push_back(&vertices[i]);
            /*std::cout << "added ";
            printVertex(vertices[i]);
            std::cout << std::endl*/
        }
    }

    return retVal;

    /*
    vector<Vertex> faceVertex;
    for (int i = 0; i < vertices.size(); ++i){
        if (vertices[i].Position.x>maxX) maxX = vertices[i].Position.x;
    }

    for (int i = 0; i < vertices.size(); ++i){
         if(abs(vertices[i].Position.x-maxX) < 0.00001){
             // we are on side with max x-values. now exclude min/max y- and z- values
             if(!floatEquals(fabs(vertices[i].Position.y), maxY) && !floatEquals(fabs(vertices[i].Position.z), maxZ))
             faceVertex.push_back(vertices [i]);
         }
    }
    return faceVertex;*/
}

vector<Vertex*> Mesh::selectVerticesMeshFaceEdge(){

    vector<Vertex*> retVal;

    for(uint i=0; i<vertices.size(); i++)
    {
        if(floatEquals(vertices[i].Position.x, maxX))
        {
            retVal.push_back(&vertices[i]);
            /*std::cout << "added ";
            printVertex(vertices[i]);
            std::cout << std::endl;*/
        }
    }

    return retVal;
}

//useful to find the index in vertices of a vertex in faceVertices
int Mesh::getIndex(Vertex* v){
    for(int i=0; i<vertices.size(); i++){
        if (vertices[i].Position==v->Position) return i;
    }

    std::cerr << "did not find vertex [" << v->Position.x << ", " << v->Position.y << ", " << v->Position.z << "]" << std::endl;
    return -1;
}

int Mesh::getIndex(int v){
    for(int i=0; i<faceVertices.size(); i++){
        if (vertices[v].Position==faceVertices[i]->Position) return i;
    }
    return -1;
}

int Mesh::getIndexTargetSurface(Vertex* v)
{
    for(int i=0; i<faceVertices.size(); i++){
        if (faceVertices[i]->Position==v->Position) return i;
    }

    //std::cerr << "did not find target-vertex [" << v->Position.x << ", " << v->Position.y << ", " << v->Position.z << "]" << std::endl;

    return -1;
}

int Mesh::getEdgeIndex(Vertex* v)
{
    for(int i=0; i<faceVerticesEdge.size(); i++){
        if (faceVerticesEdge[i]->Position==v->Position) return i;
    }

    std::cerr << "did not find vertex [" << v->Position.x << ", " << v->Position.y << ", " << v->Position.z << "]" << std::endl;
    return -1;
}

vector<int> Mesh::getNeighborsIndex(Vertex* v){
    int index = getIndex(v);
    int indexTargetSurface = getIndexTargetSurface(v);
    return getNeighborsIndex(index, indexTargetSurface);
}

vector<int> Mesh::getNeighborsIndex(int index, int indexTargetSurface)
{
    vector<int> result;
    vector<uint> faces = adjacentFaces[index];

    for(int i=0; i<faces.size(); i++){
        for(int j=0; j<3; j++){
            if (faceVertices[indexTargetSurface]->Position != vertices[indices[faces[i]][j]].Position) {
                int ind = getIndexTargetSurface(&(vertices[indices[faces[i]][j]]));
                // -1 means not found -> vertex on edge
                if(ind != -1)
                {
                    if(std::find(result.begin(), result.end(), ind) != result.end()) {
                        /* result contains ind */
                    } else {
                        /* result does not contain ind */
                        result.push_back(ind);
                    }
                }
            }
         }
    }
    return result;
}


void Mesh::createEdgeIndices()
{

    edgeIndices.clear();

    for (uint i=0; i<indices.size(); i++)
    {
        bool valid = true;

        for (uint j=0; j<3; j++)
        {
            Vertex* v = &vertices[indices[i][j]];

            if(!floatEquals(v->Position.x, maxX))
            {
                valid = false;
                break;
            }
        }

        if(valid)
        {
            edgeIndices.push_back(indices[i]);
        }
    }

    // update indices .. (currently pointing to old values)

    for(uint i=0; i<edgeIndices.size(); i++)
    {
        for (uint j=0; j<3; j++)
        {
            edgeIndices[i][j] = getEdgeIndex(&vertices[edgeIndices[i][j]]);
        }
    }
}

vector<int> Mesh::insertSorted(vector<int> vec, int in, Vertex* v2){
    vector<int> vec2=vec;
    bool under=true;
    for(int i=0; i<(vec.size()+1); i++){
        if(under){
            if(glm::distance(faceVertices[vec[i]]->Position,v2->Position)<glm::distance(faceVertices[in]->Position,v2->Position))
                vec2.push_back(vec[i]);
            else {
                vec2.push_back(in);
                under=false;
            }
        }
        else vec2.push_back(vec[i-1]);
    }
    return vec2;
}


//return the index of the 6 closest neighboors
vector<int> Mesh::getClosestNeighbors(Vertex* v){
    //vector<Vertex*> result;
    vector<int> result;
    result.push_back(0);
    glm::vec3 pivot=faceVertices[0]->Position;
    int n;

    for (int i=0; i<faceVertices.size(); i++){
        if(glm::distance(faceVertices[i]->Position,v->Position)<glm::distance(pivot,v->Position)){
            pivot=faceVertices[i]->Position;
            if (result.size()<6){
                result=insertSorted(result,i,v);
            }
            else {
                result.pop_back();
                result=insertSorted(result,i,v);
            }
        }
    }
    return result;
}

