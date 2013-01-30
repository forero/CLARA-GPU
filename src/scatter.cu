#include <cuda.h>
#include <curand_kernel.h>
#include <stdio.h>


//#include "scatter.h"

__global__ void randomStep(float *a, int N)
{

  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) a[idx] = 3.14159;

}

extern "C" void scatter_x(float *x, int min_id, int max_id){
  float *x_aux;
  float *x_aux_d;
  int n_aux;
  int blockSize;
  int nBlocks;
  int i;

  n_aux = max_id - min_id;

  if(!(x_aux = (float *)malloc(n_aux * sizeof(float)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }
  
  for(i=0;i<n_aux;i++){
    x_aux[i] = x[min_id+i];
  }

  // allocate memory on device
  cudaMalloc((void **) &x_aux_d, n_aux * sizeof(float));

  // copy data from host to device
  cudaMemcpy(x_aux_d, x_aux, sizeof(float) * n_aux, cudaMemcpyHostToDevice);

  blockSize = 32; // This is the number of threads inside a block
  nBlocks = n_aux/blockSize + (n_aux%blockSize == 0?0:1); // This is the number of blocks

  // Make the random step
  randomStep <<< nBlocks, blockSize >>> (x_aux_d, n_aux);

  // copy data from device to host
  cudaMemcpy(x_aux, x_aux_d, sizeof(float) * n_aux, cudaMemcpyDeviceToHost);
}


extern "C" void TransportPhotons(float *x, int n_photons){
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
    scatter_x(x, i*pack_size, (i+1)*pack_size);
  }
}
