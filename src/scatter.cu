#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <string.h>

#include "io.h"
#include "struct.h"
/*this makes the compiler happy*/
#include "vector.cu"


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
  FLOAT random number generator
*/
__device__ void getPoint(FLOAT *x, curandState *state)
{
  *x = curand_normal_double(state);
}


__device__ void RND_spherical(FLOAT *vec, curandState *state)
/*vector randomly distributed over the sphere*/
{
    FLOAT theta, phi;
    theta = acos(2.0*(curand_uniform(state) -0.5));
    phi = 2.0*PI*curand_uniform(state);
    vec[0] = sin(theta)*cos(phi);
    vec[1] = sin(theta)*sin(phi);
    vec[2] = cos(theta);
}

__device__ void RND_lyman_parallel_vel(FLOAT *u_parallel, FLOAT x, FLOAT a, curandState *state, int *status) 
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
    FLOAT tmp0, tmp1, tmp2;        
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

__device__ void RND_lyman_perp_vel(FLOAT *u_1, FLOAT *u_2, curandState *state)
/* 
   Genereates magnitudes for the atom's perpendicular velocities using
   Box&Muller method. 
*/
{
  FLOAT tmp1, tmp2;
  FLOAT vel_1, vel_2;
  tmp1 = curand_uniform(state);
  tmp2 = curand_uniform(state);
  vel_1 = sqrt(-log(tmp1))*cos(2.0*PI*tmp2);
  vel_2 = sqrt(-log(tmp1))*sin(2.0*PI*tmp2);
  *u_1 = vel_1;
  *u_2 = vel_2;  
  return;
}


__device__ void RND_pair(FLOAT *r_1, FLOAT *r_2, curandState *state){
  /*generates a pair of random numbers with norm less than one*/
  int finished;
  FLOAT rand_1, rand_2;

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


__device__ void RND_lyman_atom(FLOAT *Vel, FLOAT *DirPhoton, FLOAT *DirOutPhoton, FLOAT x, FLOAT a ,curandState *state, int *status)
/* Obtains a random velocity for the hydrogen atom
   the velocity is in units of the termal velocity.*/
{
    int i;
    FLOAT LocalVel[3];
    FLOAT x_axis[3];
    FLOAT y_axis[3];
    FLOAT z_axis[3];
    FLOAT rand_axis[3];
    FLOAT R_1, R_2, R_3, T, mu, iso;
    FLOAT  x_corewing;
    FLOAT exponent;

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


/*This makes the compiler happy*/
#include "propagate.cu"
#include "test.cu"

/*
  This is CLARA's core.
  This routine computes all the scatterings until the photon escapes.
*/
__global__ void scatterStep(FLOAT *x, FLOAT *p, FLOAT *k, int * n_scatter, int *status_ID, int N, curandState *const rngStates, setup *S)
{
  /*physical variables that define the photon state*/
  FLOAT pos[3];
  FLOAT dir[3];
  FLOAT r_travel, x_photon, x_comoving;
  FLOAT last_x;

  /*physical variables that define the medium*/
  FLOAT a, nu_doppler, n_HI, v_thermal, BulkVel[3], temperature;  

  /*status variables*/
  int photon_status;
  int program_status;
  int n_iter=0;
  int i;

  /*memory thread on the GPU*/
  int id = blockIdx.x*blockDim.x + threadIdx.x;
  int idx = blockIdx.x*blockDim.x + (threadIdx.x*3  + 0); 
  int idy = blockIdx.x*blockDim.x + (threadIdx.x*3  + 1);
  int idz = blockIdx.x*blockDim.x + (threadIdx.x*3  + 2);


  FLOAT f=0.0;
  FLOAT px=0.0;
  FLOAT py=0.0;
  FLOAT pz=0.0;
  FLOAT norm;


  /*Initializes the random number generator*/
  curandState localState = rngStates[id];

  /*Make the initialization for the photon*/
  pos[0] = p[idx];
  pos[1] = p[idy];
  pos[2] = p[idz];
  dir[0] = k[idx];
  dir[1] = k[idy];
  dir[2] = k[idz];
  x_photon = x[id];
  photon_status = status_ID[id]; 
  
  if (id < N){    
    while(PropagateIsInside(pos[0], pos[1], pos[2], S) && (photon_status==ACTIVE) && (n_iter<MAX_ITER)){

      /* get the temperature at this point*/
      PropagateGetTemperature(&temperature, &(pos[0]), S);
      
      /* get the number density at this point*/
      PropagateGetNumberDensity(&n_HI, pos, S);
      
      /*get the bulk velocity of the fluid at this point*/
      PropagateGetBulkVel(BulkVel, pos, S);
      
      /*Get the thermal velocity and doppler broadening*/
      nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
      a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
      v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
      
      /*change the value of the frequency to one comoving with the fluid*/
      PropagateLorentzFreqChange(&x_photon, dir, BulkVel, v_thermal, -1); 
            
      /*change the direction of the photon to the fluid frame*/
      PropagateLorentzDirChange(&(dir[0]), BulkVel, -1);
                  
      /*--------------------------------------------------------------------------*/
      /*Change the frequency and the Propagation direction, find the displacement*/	
      photon_status = PropagateStep(&x_photon, &(dir[0]), &r_travel, &a, &localState, S);	    	
      /*--------------------------------------------------------------------------*/
            
      /*Change the new direction to the lab frame value*/
      PropagateLorentzDirChange(&(dir[0]), BulkVel, 1);
      
      /*Change the frequency comoving to the lab frame value*/
      PropagateLorentzFreqChange(&x_photon, dir, BulkVel, v_thermal, 1); 
      
      /*Update the position*/
      for(i=0;i<3;i++){
	pos[i] += r_travel*dir[i];	    
      }
      
      n_iter++;
    }
    __syncthreads();
  }

  /*update the photon status*/
  if(photon_status==ACTIVE){
    photon_status = OUT_OF_BOX;
  }
  if(n_iter>= MAX_ITER){
    photon_status = SATURATED_ITERATIONS;
  }


  /*update the values*/
  p[idx] = pos[0];
  p[idy] = pos[1];
  p[idz] = pos[2];
  k[idx] = dir[0];
  k[idy] = dir[1];
  k[idz] = dir[2];
  x[id] = x_photon;  
  status_ID[id] = photon_status;
  n_scatter[id] = n_iter;

  __syncthreads();
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
extern "C" void scatter_bunch(FLOAT *x, FLOAT *p, FLOAT *k, int *n_scatter, int *status_ID, int min_id, int max_id){
  FLOAT *x_aux;
  FLOAT *p_aux;
  FLOAT *k_aux;
  int *n_scatter_aux;
  int *status_ID_aux;
  FLOAT *x_aux_d;
  FLOAT *p_aux_d;
  FLOAT *k_aux_d;
  int *n_scatter_aux_d;
  int *status_ID_aux_d;

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
  cudaError_t cudaResult = cudaSuccess;

  cudaSetDevice(0);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);

 char msg[256];
        sprintf(msg, "  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
                (float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
        printf("%s", msg);

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
  if(!(x_aux = (FLOAT *)malloc(n_aux * sizeof(FLOAT)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }

  if(!(p_aux = (FLOAT *)malloc(3 * n_aux * sizeof(FLOAT)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }

  if(!(k_aux = (FLOAT *)malloc(3 * n_aux * sizeof(FLOAT)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }

  if(!(n_scatter_aux = (int *)malloc(n_aux * sizeof(int)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }

  if(!(status_ID_aux = (int *)malloc(n_aux * sizeof(int)))){
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
    n_scatter_aux[i] = n_scatter[min_id + i];
    status_ID_aux[i] = status_ID[min_id + i];
  }


  // allocate memory on device
  cudaMalloc((void **) &x_aux_d, n_aux * sizeof(FLOAT));
  cudaMalloc((void **) &p_aux_d, 3 * n_aux * sizeof(FLOAT));
  cudaMalloc((void **) &k_aux_d, 3 * n_aux * sizeof(FLOAT));
  cudaMalloc((void **) &n_scatter_aux_d, n_aux * sizeof(int));
  cudaMalloc((void **) &status_ID_aux_d, n_aux * sizeof(int));

  // copy data from host to device
  cudaMemcpy(x_aux_d, x_aux, sizeof(FLOAT) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(p_aux_d, p_aux, 3 * sizeof(FLOAT) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(k_aux_d, k_aux, 3 * sizeof(FLOAT) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(n_scatter_aux_d, n_scatter_aux, sizeof(int) * n_aux, cudaMemcpyHostToDevice);
  cudaMemcpy(status_ID_aux_d, status_ID_aux, sizeof(int) * n_aux, cudaMemcpyHostToDevice);



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
  scatterStep <<< nBlocks, blockSize >>> (x_aux_d, p_aux_d, k_aux_d, n_scatter_aux_d, status_ID_aux_d, n_aux, d_rngStates, S);

  // copy data from device to host
  cudaMemcpy(x_aux, x_aux_d, sizeof(FLOAT) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(p_aux, p_aux_d, 3 * sizeof(FLOAT) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(k_aux, k_aux_d, 3 * sizeof(FLOAT) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(n_scatter_aux, n_scatter_aux_d, sizeof(int) * n_aux, cudaMemcpyDeviceToHost);
  cudaMemcpy(status_ID_aux, status_ID_aux_d, sizeof(int) * n_aux, cudaMemcpyDeviceToHost);

  printf("%f\n", All.Tau);

  for(i=0;i<n_aux;i++){
    x[min_id+i] = x_aux[i];
    for(l=0;l<3;l++){
      p[3*(min_id+i)+l] = p_aux[3*i + l]; 
      k[3*(min_id+i)+l] = k_aux[3*i + l]; 
    }    
    n_scatter[min_id + i] = n_scatter_aux[i];
    status_ID[min_id + i] = status_ID_aux[i];
  }
}

