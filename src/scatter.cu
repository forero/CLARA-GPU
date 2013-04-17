#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "struct.h"


/*
  RNG kernel initialization
*/
__global__ void initRNG(curandState *const rngStates,
                        const unsigned int seed)
{
    // Determine thread ID
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Initialise the RNG
    curand_init(seed, tid, 0, &rngStates[tid]);
}

/*
  float random number generator
*/
__device__ void getPoint(float *x, curandState *state)
{
  *x = curand_normal_double(state);
}


/* 
   Generates the parallel velocity to the photon propagation.

   First generates a random number \theta between -pi/2 and pi/2, then 
   it generates u_parallel through u_parallel = a \tan\theta + x
   then this value of u_parallel is kept if another random number 
   between [0,1] is smaller than \exp^{-u_parallel^2}. 
   Finally, the value of u_parallel is multiplied by the sign of x.       
*/
__device__ void RND_lyman_parallel_vel(float *u_parallel, float x, float a, curandState *state, int *status) 
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

/*
  Geometrical test. It decides whether the photon is still inside the
  propagation volume.
*/
__device__ int PropagateIsInside(float pos_x, float pos_y, float pos_z, setup *S)
{
    int is_in;
    double radius;
    is_in = 1;

    /*Infinite Slab*/
    if(S->NeufeldSlab){
	if(fabs(pos_z)<S->Tau){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    /*Cube*/
    if(S->NeufeldCube){
	if(fabs(pos_x)<S->Tau && 
	   fabs(pos_y)<S->Tau &&
	   fabs(pos_z)<S->Tau){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    /*Sphere*/
    if(S->ExpandingSphere||S->RotatingSphere){
      radius = pos_x*pos_x + pos_y*pos_y + pos_z*pos_z;
      radius = sqrt(radius);
      if(radius< S->Tau){
	is_in = 1;
      }else{
	is_in = 0 ;
      }
    }
    
    return is_in;
}


/*
  This is CLARA's core.
  This routine computes all the scatterings until the photon escapes.
*/
__global__ void scatterStep(float *x, float *p, float *k, int N, curandState *const rngStates, setup *S)
{
  /*physical variables that define the medium*/
  float HIInteractionProb, dust_sigma;
  float nu_doppler, v_thermal, temperature;

  /*random variables that will define dust or hydrogen interaction*/
  float rand_interaction, rand_absorption;
  
  /*physical variables that defie the photon state*/
  float v_parallel, g_recoil, lya_sigma, tau, x_in, x_out, r, a_in;
  float k_out_photon[3];
  float u_atom[3];

  /*memory thread on the GPU*/
  int id = blockIdx.x*blockDim.x + threadIdx.x;
  int idx = blockIdx.x*blockDim.x + (threadIdx.x*3  + 0); 
  int idy = blockIdx.x*blockDim.x + (threadIdx.x*3  + 1);
  int idz = blockIdx.x*blockDim.x + (threadIdx.x*3  + 2);


  float f=0.0;
  float px=0.0;
  float py=0.0;
  float pz=0.0;
  float norm;
  int status;


  /*Initializes the random number generator*/
  curandState localState = rngStates[id];


  if (id < N){    
    while(PropagateIsInside(p[idx], p[idy], p[idz], S)){
      RND_lyman_parallel_vel(&f, x[id], 0.01, &localState, &status); 

      getPoint(&px, &localState);
      getPoint(&py, &localState);
      getPoint(&pz, &localState);
            
      p[idx] = p[idx] + px;
      p[idy] = p[idy] + py;
      p[idz] = p[idz] + pz;
      x[id] = x[id] + f/1E5;          
      __syncthreads();
    }
  }
  
}

/*
  This is the main driver for CLARA.
  Initializes the memory on the device for:
  - Photons positions.
  - Photons direction of propagation.
  - Photons frequency.
  - Global setup (densities, velocities, temperatures)

  This is also the place where the main GPU characteristics have to be setup:
  - Number of blocks
  - Number of threads
*/
extern "C" void scatter_bunch(float *x, float *p, float *k, int min_id, int max_id){
  float *x_aux;
  float *p_aux;
  float *k_aux;
  float *x_aux_d;
  float *p_aux_d;
  float *k_aux_d;

  int n_aux;
  int blockSize;
  int nBlocks;
  int i,l;
  setup *S;
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

  //allocate and copy the global setup structure 
  cudaMalloc((void**) &S, sizeof(setup));
  cudaMemcpy(S, &All, sizeof(setup), cudaMemcpyHostToDevice);

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

  if(!(k_aux = (float *)malloc(3 * n_aux * sizeof(float)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }

  //fill the auxiliary variables
  for(i=0;i<n_aux;i++){
    x_aux[i] = x[min_id+i];
    for(l=0;l<3;l++){
      p_aux[3*i + l] = p[3*(min_id+i)+l];
      k_aux[3*i + l] = k[3*(min_id+i)+l];
    }
  }


  // allocate memory on device
  cudaMalloc((void **) &x_aux_d, n_aux * sizeof(float));
  cudaMalloc((void **) &p_aux_d, 3 * n_aux * sizeof(float));
  cudaMalloc((void **) &k_aux_d, 3 * n_aux * sizeof(float));

  // copy data from host to device
  cudaMemcpy(x_aux_d, x_aux, sizeof(float) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(p_aux_d, p_aux, 3 * sizeof(float) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(k_aux_d, k_aux, 3 * sizeof(float) * n_aux, cudaMemcpyHostToDevice);

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
  scatterStep <<< nBlocks, blockSize >>> (x_aux_d, p_aux_d, k_aux_d, n_aux, d_rngStates, S);

  // copy data from device to host
  cudaMemcpy(x_aux, x_aux_d, sizeof(float) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(p_aux, p_aux_d, 3 * sizeof(float) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(k_aux, k_aux_d, 3 * sizeof(float) * n_aux, cudaMemcpyDeviceToHost);

  printf("%f\n", All.Tau);

  for(i=0;i<n_aux;i++){
    x[min_id+i] = x_aux[i];
    for(l=0;l<3;l++){
      p[3*(min_id+i)+l] = p_aux[3*i + l]; 
      k[3*(min_id+i)+l] = k_aux[3*i + l]; 
    }    
  }
}
