#include "utils.h"

int outTriplet(vector<int> vec, int begin, int end){
    for (int i=begin; i<=end; i++){
         for (int j=begin; j<=end; j++)
            if (i!=vec[j]) return i;

    }
}

bool floatEquals(float val1, float val2)
{
    return (fabs(val1 - val2) < 0.00001);
}

float fbar(float x, float dth){
    float f= fmax(0, -std::log((1-x)+dth));
    return f;
}

glm::vec3 proj(glm::vec3 xs, glm::vec3 di, glm::vec3 pos){
    glm::vec3 result;
    result.x=pos.x;
    result.y=xs.y;
    result.z=xs.z;
    return result;
}

float* matrixProduct(array* L, float* X){
    float* R = new float[NORMALS];
    for(int i=0; i<NORMALS; i++){
        for(int j=0; j<NORMALS; j++)
            R[i]+=L[i][j]*X[j];
    }
    return R;
}

void printMatrix(array* X){
    for(int i=0; i<NORMALS; i++){
        std::cout<<"| "<<std::flush;
        for(int j=0; j<NORMALS; j++)
            std::cout<<" "<<X[i][j]<<" "<<std::flush;
        std::cout<<" |"<<std::endl;
    }

}

bool isEdge(glm::vec3 v){
    bool edge = false;
    if(floatEquals(fabs(v[1]), MAX_Y) || floatEquals(fabs(v[2]), MAX_Z)) edge=true;
    return edge;

}

bool compareVectors(vector<int> v1, vector<int> v2){
    bool same = true;
    if(v1.size()!=v2.size()){
        std::cout<<"v1 size: "<<v1.size()<<std::endl;
        std::cout<<"v2 size: "<<v2.size()<<std::endl;
        return false;
    }
    for(uint i = 0; i<v1.size(); i++){
        if(v1[i]!=v2[i]){
            std::cout<<"different vectors"<<std::endl;
            return false;
        }
    }
}
