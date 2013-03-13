#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "struct.h"


// RNG init kernel
__global__ void initRNG(curandState *const rngStates,
                        const unsigned int seed)
{
    // Determine thread ID
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Initialise the RNG
    curand_init(seed, tid, 0, &rngStates[tid]);
}

__device__ void getPoint(float *x, curandState *state)
{
  *x = curand_normal_double(state);
}


__device__ void RND_lyman_parallel_vel(float *u_parallel, float x, float a, curandState *state, int *status) 
/* 
   it generates a random number \theta between -pi/2 and pi/2, then 
   it generates u_parallel through u_parallel = a \tan\theta + x
   then this value of u_parallel is kept if another random number 
   between [0,1] is smaller than \exp^{-u_parallel^2}. 
   At the end the value of u_parallel is multiplied by the sign of x.       
*/
{
    int finished = 0;
    float tmp0, tmp1, tmp2;        
    int counter = 0;
    
    while (finished == 0) {
      tmp0 = (curand_uniform(state) - 0.5)*PI ;
      tmp1 = (a*tan(tmp0)) + fabs(x); 
      tmp2  = curand_uniform(state);
      if(tmp2 <= (exp(-(tmp1*tmp1)))) finished = 1;		
      counter++;		
      if(counter > MAX_VEL_ITER) {
	finished = 1;		    
	*status = EXCEEDED_ITERATIONS;
      }	
    }

    if(x > 0.0){
      *u_parallel = tmp1;
    }else{    
      *u_parallel = -tmp1;
    }
}


__global__ void scatterStep(float *x, float *p, int N, curandState *const rngStates)
{

  int id = blockIdx.x*blockDim.x + threadIdx.x;
  int idx = blockIdx.x*blockDim.x + (threadIdx.x*3  + 0);
  int idy = blockIdx.x*blockDim.x + (threadIdx.x*3  + 1);
  int idz = blockIdx.x*blockDim.x + (threadIdx.x*3  + 2);

  int status;

  // Initialise the RNG
  curandState localState = rngStates[id];
  float f=0.0;
  float px=0.0;
  float py=0.0;
  float pz=0.0;
  float norm;

  if (id < N){    
    while((fabs(p[idx])<5.0) || (fabs(p[idy])<10.0) || (fabs(p[idz])<15.0)){
      RND_lyman_parallel_vel(&f, x[id], 0.01, &localState, &status); 

      getPoint(&px, &localState);
      getPoint(&py, &localState);
      getPoint(&pz, &localState);
      
      p[idx] = p[idx] + px;
      p[idy] = p[idy] + py;
      p[idz] = p[idz] + pz;
      x[id] = x[id] + f/1000.0;          
      __syncthreads();
    }
  }
  
}

extern "C" void scatter_bunch(float *x, float *p, int min_id, int max_id){
  float *x_aux;
  float *p_aux;
  float *x_aux_d;
  float *p_aux_d;
  int n_aux;
  int blockSize;
  int nBlocks;
  int i,k;

  unsigned int m_seed;
  unsigned int m_numSims;
  unsigned int m_device;
  unsigned int m_threadBlockSize;  
  struct cudaDeviceProp     deviceProperties;
  struct cudaFuncAttributes funcAttributes;
  cudaError_t cudaResult = cudaSuccess;

  cudaResult = cudaGetLastError();
  if (cudaResult!=cudaSuccess){
    printf( "Error after device!\n" );
    printf("CUDA error: %s\n", cudaGetErrorString(cudaResult));
  }  

  
  m_device = 0;
  // Get device properties
  cudaResult = cudaGetDeviceProperties(&deviceProperties, m_device);
  fprintf(stdout, "Device Properties:\n Multiproc count %d\n", deviceProperties.multiProcessorCount );

  //allocate auxiliary variables  
  n_aux = max_id - min_id;
  if(!(x_aux = (float *)malloc(n_aux * sizeof(float)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }

  if(!(p_aux = (float *)malloc(3 * n_aux * sizeof(float)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }
  
  for(i=0;i<n_aux;i++){
    x_aux[i] = x[min_id+i];
    for(k=0;k<3;k++){
      p_aux[3*i + k] = p[3*(min_id+i)+k];
    }
  }

  // allocate memory on device
  cudaMalloc((void **) &x_aux_d, n_aux * sizeof(float));
  cudaMalloc((void **) &p_aux_d, 3 * n_aux * sizeof(float));

  // copy data from host to device
  cudaMemcpy(x_aux_d, x_aux, sizeof(float) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(p_aux_d, p_aux, 3 * sizeof(float) * n_aux, cudaMemcpyHostToDevice);

  blockSize = 512; // This is the number of threads inside a block
  nBlocks = (3*n_aux)/blockSize + (n_aux%blockSize == 0?0:1); // This is the number of blocks
  fprintf(stdout, "nBlocks %d\n", nBlocks);

  //allocate memory for RNG states
  curandState *d_rngStates = 0;
  cudaResult = cudaMalloc((void **)&d_rngStates, blockSize * nBlocks * sizeof(curandState));

  // Initialise RNG
  m_seed = min_id;
  initRNG<<<nBlocks, blockSize>>>(d_rngStates, m_seed);

  // Make the random step until all the photons escape
  scatterStep <<< nBlocks, blockSize >>> (x_aux_d, p_aux_d, n_aux, d_rngStates);

  // copy data from device to host
  cudaMemcpy(x_aux, x_aux_d, sizeof(float) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(p_aux, p_aux_d, 3 * sizeof(float) * n_aux, cudaMemcpyDeviceToHost);

  for(i=0;i<n_aux;i++){
    x[min_id+i] = x_aux[i];
    for(k=0;k<3;k++){
      p[3*(min_id+i)+k] = p_aux[3*i + k]; 
    }    
  }
}
