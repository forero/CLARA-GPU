#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "struct.h"
#include "vector.h"

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


__device__ void RND_spherical(float *vec, curandState *state)
/*vector randomly distributed over the sphere*/
{
    float theta, phi;
    theta = acos(2.0*(curand_uniform(state) -0.5));
    phi = 2.0*PI*curand_uniform(state);
    vec[0] = sin(theta)*cos(phi);
    vec[1] = sin(theta)*sin(phi);
    vec[2] = cos(theta);
}

__device__ void RND_lyman_parallel_vel(float *u_parallel, float x, float a, curandState *state, int *status) 
/* 
   Generates the parallel velocity to the photon propagation.
   First generates a random number \theta between -pi/2 and pi/2, then 
   it generates u_parallel through u_parallel = a \tan\theta + x
   then this value of u_parallel is kept if another random number 
   between [0,1] is smaller than \exp^{-u_parallel^2}. 
   Finally, the value of u_parallel is multiplied by the sign of x.       
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

__device__ void RND_lyman_perp_vel(float *u_1, float *u_2, curandState *state)
/* 
   Genereates magnitudes for the atom's perpendicular velocities using
   Box&Muller method. 
*/
{
  float tmp1, tmp2;
  float vel_1, vel_2;
  tmp1 = curand_uniform(state);
  tmp2 = curand_uniform(state);
  vel_1 = sqrt(-log(tmp1))*cos(2.0*PI*tmp2);
  vel_2 = sqrt(-log(tmp1))*sin(2.0*PI*tmp2);
  *u_1 = vel_1;
  *u_2 = vel_2;  
  return;
}


__device__ void RND_pair(float *r_1, float *r_2, curandState *state){
  /*generates a pair of random numbers with norm less than one*/
  int finished;
  float rand_1, rand_2;

  finished = 0;
  while (finished == 0) {
    rand_1 = 2.0*(curand_uniform(state) - 0.5);
    rand_2 = 2.0*(curand_uniform(state) - 0.5);
    if(rand_1*rand_1 + rand_2*rand_2 < 1.0){
      finished = 1;
    }
  }    
  *r_1 = rand_1;
  *r_2 = rand_2;
}


__device__ void RND_lyman_atom(float *Vel, float *DirPhoton, float *DirOutPhoton, float x, float a ,curandState *state, int *status)
/* Obtains a random velocity for the hydrogen atom
   the velocity is in units of the termal velocity.*/
{
    int i;
    float LocalVel[3];
    float x_axis[3];
    float y_axis[3];
    float z_axis[3];
    float rand_axis[3];
    float R_1, R_2, R_3, T, mu, iso;
    float  x_corewing;
    float exponent;

    /*get first the parallel velocity*/
    RND_lyman_parallel_vel(&(LocalVel[2]), x, a, state, status);

    /*get the perpendicular velocity*/
    RND_lyman_perp_vel(&(LocalVel[0]), &(LocalVel[1]), state);

    /*get the axis in the coordinate system of the atom, where 
      the z direction is the propagation direction of the photon*/
    for(i=0;i<3;i++){
	z_axis[i] = DirPhoton[i];
    }

    /*get another random vector*/
    RND_spherical(rand_axis, state);

    /*make the cross product and get y_axis*/
    cross_product(z_axis, rand_axis, y_axis);

    /*make the cross product and get x_axis*/
    cross_product(y_axis, z_axis, x_axis);

    /*normalize the vectors*/
    normalize(x_axis);
    normalize(y_axis);
    normalize(z_axis);

    /*see if they are perpendicular*/
    rand_axis[0] = point_product(x_axis, z_axis);
    rand_axis[1] = point_product(x_axis, y_axis);
    rand_axis[2] = point_product(y_axis, z_axis);
    if(fabs(rand_axis[0]) + fabs(rand_axis[1])+ fabs(rand_axis[2])>1.0e-10){
      *status = NULL_NORM;
    }

    /*Now make the transformation into the coordinate frame of the lab*/
    for(i=0;i<3;i++){
	Vel[i] = LocalVel[0]*x_axis[i] + LocalVel[1]*y_axis[i] + LocalVel[2]*z_axis[i];
    }

    /*
      now get the outgoing direction of the photon, 
      taking advantage of the vector basis I have just generated
      Here I just take the value of dijkstra for the wing.
    */

    /*first define if it's in the core or not*/    

    x_corewing = 1.59 - 0.60*log10(a) - 0.03*log10(a)*log10(a);

    R_1 = curand_uniform(state);

    iso = curand_uniform(state);
    if(iso<(1.0/3.0)){/*isotropic*/
      mu = (2.0*R_1 - 1.0);
    }else{
      if(fabs(x)<x_corewing){			/*In the core*/
	/*now we make the decision if it's isotropic or not*/
	T = (1.0/7.0)*(14.0  - 24.0*R_1  + sqrt(245.0 - 672.0*R_1 + 576*R_1*R_1));	  
	exponent = 1.0/3.0; /*this makes cuda compiler happy*/
	mu = 1.0/(pow(T,exponent)) - pow(T, exponent);	
      }else{
	T = 2.0 - 4.0*R_1  + sqrt(5.0 -16.0*R_1 + 16*R_1*R_1);
	exponent = 1.0/3.0; /*this makes cuda compiler happy*/
	mu = 1.0/(pow(T,exponent)) - pow(T, exponent);	
      }
    }
    
    RND_pair(&R_1, &R_2, state);    
    R_3 = R_1*R_1 + R_2*R_2;
    for(i=0;i<3;i++){
	DirOutPhoton[i] = 
	    sqrt((1.0-(mu*mu))/R_3)*R_1*x_axis[i] + 
	    sqrt((1.0-(mu*mu))/R_3)*R_2*y_axis[i] + 
	    mu*z_axis[i];
    }
}


__device__ int PropagateIsInside(float pos_x, float pos_y, float pos_z, setup *S)
/*
  Geometrical test. It decides whether the photon is still inside the
  propagation volume.
*/
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
