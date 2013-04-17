#ifndef VECTOR_H
#define VECTOR_H
__device__ void cross_product(float *vec_1, float *vec_2, float *result);
__device__ float point_product(float *vec_1, float *vec_2);
__device__ void normalize(float *vec);
#endif
