#include <stdio.h>
#include <math.h>

__device__ void TestFirstScatter(float *x, setup *S, curandState *state){
    FLOAT x_in;
    FLOAT r_travel;
    FLOAT k_in_photon[3];
    FLOAT n_HI;
    int status;
    int i;
    FLOAT FileName[MAX_FILENAME_SIZE];
    FILE *out;

    int id = blockIdx.x*blockDim.x + threadIdx.x;

    x_in = x[id];
    n_HI = S->NumberDensityHI;
    RND_spherical(&(k_in_photon[0]));    
    status = PropagateStep(&x_in, &(k_in_photon), &tau_travel, a, state, S);
    x[id] = x_in;
} 
