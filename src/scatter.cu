#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>



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


__global__ void scatterStep(float *x, float *p, int N, curandState *const rngStates)
{

  int id = blockIdx.x*blockDim.x + threadIdx.x;
  int idx = blockIdx.x*blockDim.x + (threadIdx.x*3 + 0);
  int idy = blockIdx.x*blockDim.x + (threadIdx.x*3 + 1);
  int idz = blockIdx.x*blockDim.x + (threadIdx.x*3 + 2);

  // Initialise the RNG
  curandState localState = rngStates[id];
  float f=0.0;
  float px=0.0;
  float py=0.0;
  float pz=0.0;
  float norm;

  if (id < N){    
    getPoint(&f, &localState);
    getPoint(&px, &localState);
    getPoint(&py, &localState);
    getPoint(&pz, &localState);

    while(fabs(x[id])<10.0){
      p[idx] = p[idx] + px;
      p[idy] = p[idy] + py;
      p[idz] = p[idz] + pz;
      x[id] = x[id] + f;          
      getPoint(&f, &localState);
      getPoint(&px, &localState);
      getPoint(&py, &localState);
      getPoint(&pz, &localState);
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

  blockSize = 128; // This is the number of threads inside a block
  nBlocks = n_aux/blockSize + (n_aux%blockSize == 0?0:1); // This is the number of blocks

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


extern "C" void TransportPhotons(float *x, float *p, float *k, int n_photons){
  int pack_size, last_pack_size;  
  int n_packs;
  int i;

  pack_size = 400;
  
  n_packs = n_photons/pack_size;
  last_pack_size = n_photons%pack_size;
  
  fprintf(stdout, "Photons to transport %d\n", n_photons);
  fprintf(stdout, "%d packs of size %d\n", n_packs, pack_size);
  fprintf(stdout, "one last pack of size %d\n", last_pack_size);
  
  for(i=0;i<n_packs;i++){
    scatter_bunch(x, p, i*pack_size, (i+1)*pack_size);
  }  

  scatter_bunch(x, p, i*pack_size, i*pack_size + last_pack_size);  

}
