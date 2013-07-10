#include <stdio.h>
#include <math.h>
#include "struct.h"

__device__ FLOAT point_product(FLOAT *vec_1, FLOAT *vec_2){
    double point;
    int i;
    point = 0.0;
    for(i=0;i<3;i++){
        point+=vec_1[i]*vec_2[i];
    }
    return point;
}


__device__ void cross_product(FLOAT *vec_1, FLOAT *vec_2, FLOAT *result){
    result[0] = vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1];
    result[1] = vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2];
    result[2] = vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0];
    return;
}


__device__ void normalize(FLOAT *vec){
    FLOAT norm;
    int i;
    norm = 0.0;
    norm = point_product(vec, vec);
    norm = sqrt(norm);
    for(i=0;i<3;i++){
        vec[i] = vec[i]/norm;
    }
    return;
}
