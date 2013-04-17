#include "struct.h"
#ifndef VECTOR_H
#define VECTOR_H
__device__ void cross_product(FLOAT *vec_1, FLOAT *vec_2, FLOAT *result);
__device__ FLOAT point_product(FLOAT *vec_1, FLOAT *vec_2);
__device__ void normalize(FLOAT *vec);
#endif
