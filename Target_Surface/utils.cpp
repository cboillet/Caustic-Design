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
