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

