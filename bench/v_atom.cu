#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <string.h>
#include "../src/struct.h"
#include "../src/vector.cu"
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

__global__ void generate_velocites(FLOAT x, FLOAT a, curandState *states){
    // Determine thread ID
    FLOAT vel[3];
    FLOAT k_in[3];
    FLOAT k_out[3];
    int i;
    curandState state;
    int status;
    unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;

    // Initialise the RNG
    state = states[id];
    
    for(i=0;i<3;i++){
      vel[i] = 0; k_in[i]=0.0; k_out[i]=0.0;
    }
    vel[0]=1;k_in[0]=1;k_out[0]=1;
    RND_lyman_atom(&(vel[0]), &(k_in[0]), &(k_out[0]), x, a ,&state, &status);  
}



int main(int argc, char **argv){
  int blockSize, nBlocks;
  int n_points;
  curandState *d_rngStates = 0;
  cudaError_t cudaResult = cudaSuccess;
  FLOAT x, a;

  if(argc!=4){
    fprintf(stderr, "USAGE: ./a.out n_points x a \n");
    exit(1);
  }
  n_points = atoi(argv[1]);
  x = atof(argv[2]);
  a = atof(argv[3]);

  fprintf(stdout, "The number of points to generate %d, x= %f, a=%f\n", n_points, x, a);

  /*GPU blocks and threads*/
  blockSize = 512; // This is the number of threads inside a block
  nBlocks = (3 * n_points)/blockSize + (n_points%blockSize == 0?0:1); // This is the number of blocks
  fprintf(stdout, "nBlocks %d\n", nBlocks);

  /*initalize random number generator*/
  cudaResult = cudaMalloc((void **)&d_rngStates, blockSize * nBlocks * sizeof(curandState));
  initRNG<<<nBlocks, blockSize>>>(d_rngStates, 1000);

  cudaResult = cudaGetLastError();
  if (cudaResult!=cudaSuccess){
    fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(cudaResult));
    exit(1);
  }  
  generate_velocites<<<nBlocks, blockSize>>>(x, a, d_rngStates);  
  return 0;
}
