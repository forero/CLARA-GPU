#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "struct.h"
#include "ray.h"
//#include "RND_lyman.h"
#include "io.h"
#include "myrand.h"


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

extern "C" void get_lyman_perp_vel(float *u_1, float *u_2, int n_points)
/* 
   Genereates magnitudes for the atom's perpendicular velocities using
   Box&Muller method. 
*/
{

  float tmp1, tmp2;
  float vel_1, vel_2;
  int i;
  for(i=0;i<n_points;i++){
    tmp1 = drand48();
    tmp2 = drand48();
    vel_1 = sqrt(-log(tmp1))*cos(2.0*PI*tmp2);
    vel_2 = sqrt(-log(tmp1))*sin(2.0*PI*tmp2);
    u_1[i] = vel_1;
    u_2[i] = vel_2;  
  }

  return;
}

__global__ void get_lyman_parallel_vel(float *u_parallel, float *x, double a, curandState *state, int n_points) 
/* 
   Generates the parallel velocity to the photon propagation.
   First generates a random number \theta between -pi/2 and pi/2, then 
   it generates u_parallel through u_parallel = a \tan\theta + x
   then this value of u_parallel is kept if another random number 
   between [0,1] is smaller than \exp^{-u_parallel^2}. 
   Finally, the value of u_parallel is multiplied by the sign of x.       
*/
{
  int id = blockIdx.x*blockDim.x + threadIdx.x;
  int counter=0;
  int finished=0;
  float tmp0, tmp1, tmp2;
  __shared__ float my_a;
  __shared__ float my_x;

  my_a = a;
  my_x = x[id];
  
  if(id<n_points){
    while (finished == 0) {
      tmp0 = (curand_uniform(&(state[id])) - 0.5)*PI ;
      tmp1 = (my_a * tan(tmp0)) + fabs(my_x); 
      tmp2  = curand_uniform(&(state[id]));
      if(tmp2 <= (exp(-(tmp1*tmp1)))) finished = 1;		
      counter++;		
      if(counter > MAX_VEL_ITER) {
	finished = 1;		    
      }	
    }
    
    if(x[id] > 0.0){
      u_parallel[id] = tmp1;
    }else{    
      u_parallel[id] = -tmp1;
    }
  }
  //  __syncthreads();
}


extern "C" double point_product(double *vec_1, double *vec_2){
    double point;
    int i;
    point = 0.0;
    for(i=0;i<3;i++){
        point+=vec_1[i]*vec_2[i];
    }
    return point;
}


extern "C" void cross_product(double *vec_1, double *vec_2, double *result){
    result[0] = vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1];
    result[1] = vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2];
    result[2] = vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0];
}



extern "C" void normalize(double *vec){
    double norm;
    int i;
    norm = 0.0;
    norm = point_product(vec, vec);
    norm = sqrt(norm);
    for(i=0;i<3;i++){
        vec[i] = vec[i]/norm;
    }
}



extern "C" double LyaTau(void)
/*generates a value of optical depth*/
{
    double tau;
    double tmp;
    tmp = drand48();

    tau = -log(tmp);
    return tau;
}

extern "C" double LyaH(double x, double a)
{
    double z, P, q, H;
    z = 0.0; P=0.0; q=0.0; H=0.0;

    z = (x*x - 0.855)/(x*x + 3.42);
    
    P = 5.674*(z*z*z*z) - 9.207*(z*z*z)  + 4.421*(z*z) + 0.1117*z;

    if(z<=0)
    {
	q = 0.0;
    }
    else
    {
	q = (1.0 + (21.0/(x*x)))*(a/(PI*(x*x + 1.0)))*P;
    }

    H = q*sqrt(PI);
    H = H + exp(-(x*x));
    return H;
}

extern "C" double LyaCrossSection(double x, double a){
    double sigma;
    double nu_doppler;


    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    sigma = 0.4162*sqrt(PI)*(CHARGEELECTRON*CHARGEELECTRON)/(ELECTRONMASS*C_LIGHT*nu_doppler);
    sigma = sigma*LyaH(x,a);
        
    return sigma;
}


#include "RND_lyman.c"

extern "C" int PropagateIsInside(double PosX, double PosY, double PosZ)
/*depending on the geometrical consideration of the problem at hand,
  decides if the photon is still inside*/
{
    int is_in;
    double radius;
    is_in = 1;

    //    fprintf(stdout, "x, y, z i _d %f %f %f %f\n", PosX, PosY, PosZ, All.SlabLength);

    if(All.SimulationCube){
	if(fabs(PosX)<All.SlabLength && 
	   fabs(PosY)<All.SlabLength &&
	   PosZ < (All.SlabLength*2.0) && PosZ>0.0){
	  is_in = 1;
	}else{
	  is_in = 0 ;
	}
    }

    if(All.NeufeldSlab){
	if(fabs(PosX)<All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.NeufeldCube){
	if(fabs(PosX)<All.SlabLength && 
	   fabs(PosY)<All.SlabLength &&
	   fabs(PosZ)<All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }
    
    if(All.ExpandingSphere){
	radius = PosX*PosX + PosY*PosY + PosZ*PosZ;
	radius = sqrt(radius);
	if(radius< All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }
    
    return is_in;
}

void PropagateGetBulkVel(double *BulkVel, double *Pos){
    int i;

    for(i=0;i<3;i++){
	BulkVel[i] = 0.0;
    }

    if(All.ExpandingSphere){
	for(i=0;i<3;i++){
	    BulkVel[i] = (Pos[i]/All.SlabLength)*(All.VmaxSphere*1.0e5); /*in cm/s*/
	}
    }
}

void PropagateGetTemperature(double *Temp, double *Pos){
    *Temp = All.Temperature;
}

void PropagateGetNumberDensity(double *n_HI, double *Pos){
    *n_HI = All.NumberDensityHI;
}


void PropagateLorentzDirChange(double *Dir, double *BulkVel, int sign){
    int i;
    for(i=0;i<3;i++){
	Dir[i] = Dir[i]*(1.0 + sign*(BulkVel[i]/C_LIGHT)); 
    }
    
    /*renormalize, I just want to know the direction change, bust still be normalized*/
    normalize(Dir);
}

void PropagateLorentzFreqChange(float *x, double *Dir, 
			    double *BulkVel, double ThermalVel, int sign){
    float lorentz_factor;
    lorentz_factor = (BulkVel[0]*Dir[0] + BulkVel[1]*Dir[1] + BulkVel[2]*Dir[2])/ThermalVel;
    
    *x  = *x + sign*lorentz_factor;
}

int count_active(int *status, int n_points){
  int i, n_active;
  n_active = 0;
  for(i=0;i<n_points;i++){
    if(status[i]==ACTIVE){
      n_active++;
    }
  }
  return n_active;
}

void PhotonInitialize(double *PosX, double *PosY, double *PosZ, double *DirX, double *DirY, double *DirZ){
    int i;
    double theta, phi;


    RND_spherical(DirX, DirY, DirZ);
    *PosX = 0.0;
    *PosY = 0.0;
    *PosZ = 0.0;
    
    if(All.HomogeneousInit){
      if(All.NeufeldCube){
	*PosX = 2.0*(drand48()-0.5)*All.SlabLength;
	*PosY = 2.0*(drand48()-0.5)*All.SlabLength;
	*PosZ = 2.0*(drand48()-0.5)*All.SlabLength;
      }
      
      if(All.ExpandingSphere){
	*PosX = 2.0*(drand48()-0.5)*All.SlabLength;
	*PosY = 2.0*(drand48()-0.5)*All.SlabLength;
	*PosZ = 2.0*(drand48()-0.5)*All.SlabLength;
	do{
	  *PosX = 2.0*(drand48()-0.5)*All.SlabLength;
	  *PosY = 2.0*(drand48()-0.5)*All.SlabLength;
	  *PosZ = 2.0*(drand48()-0.5)*All.SlabLength;
	}while(PropagateIsInside(*PosX, *PosY, *PosZ)==0);	    
      }	  
      
      if(All.NeufeldSlab){
	*PosX = 2.0*(drand48()-0.5)*All.SlabLength;
	*PosY = 0.0;
	*PosZ = 0.0;
      }
    }    

    if(All.SimulationCube){
      theta = acos(drand48());
      phi = 2.0*PI*drand48();
      *DirX = sin(theta)*cos(phi);
      *DirY = sin(theta)*sin(phi);
      *DirZ = cos(theta);
      *PosX = 2.0*(drand48()-0.5)*All.SlabLength;
      *PosY = 2.0*(drand48()-0.5)*All.SlabLength;
      *PosZ = All.SlabLength/All.Tau;
    }

}



int PropagateStep(float *v_parallel, float *v_perp_1, float *v_perp_2, float *x_in, float *x_out, double *k_in_x, double *k_in_y, double *k_in_z, double *r_travel, int *status, double a, double n_HI, int n_points) 
/* Calculates: 
   1. new frequency 
   2. new propagation direction, 
   3. distance to travel to next interaction

   Note:
   At this point the frequency must be already the comoving one (not lab frame!)
   If the photon is absorbed by dust, then the status is set to -1;
*/
{
    double lya_sigma;
    double tau;
    double r;
    double N;
    double a_in;
    int i;
    double HIInteractionProb; /*probability of interacting with HI*/
    double dust_sigma;
    double rand_interaction;
    double rand_absorption;
    double nu_doppler;
    double g_recoil;
    double temperature;
    int i_photon;


    /*Initalizations*/
    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    temperature = (nu_doppler/CONSTANT_NU_DOPPLER)*(nu_doppler/CONSTANT_NU_DOPPLER)*10000.0;
    g_recoil = PLANCK*nu_doppler/(2.0*BOLTZMANN*temperature);

    /* generates the atom velocity the comoving frame*/
    RND_lyman_atom(v_parallel, v_perp_1, v_perp_2, k_in_x, k_in_y, k_in_z, x_out, status, a, g_recoil, n_points);
    
    for(i_photon=0;i_photon<n_points;i_photon++){
      if(status[i_photon]==ACTIVE){

	HIInteractionProb = 1.0;
	rand_absorption = 0.0;
	
	/*get a random number to decide if the thing interact with HI*/
	rand_interaction  = drand48();      
	N    = n_HI;

	/* basic status at the interacion point*/
	lya_sigma  = LyaCrossSection(x_in[i_photon], a);
	tau        = LyaTau();
	r          = tau/(lya_sigma*N); /*the distance is in cm*/

	if(All.UseDust){
	  dust_sigma = PI*All.GrainSize*All.GrainSize;
	  All.NumberDensityDust = All.TauDust/(All.SlabLength*dust_sigma);
	  /*probability of interacting with HI*/
	  HIInteractionProb = 
	    All.NumberDensityHI*lya_sigma/(All.NumberDensityHI*lya_sigma + All.NumberDensityDust*dust_sigma);
	  /*	fprintf(stdout, "Interaction prob %e \n", HIInteractionProb);*/
	  /*get a random number to decide if the thing will be absorbed*/
	  rand_absorption = drand48();	  
	  /* The traveled distance is modifyed by the presence of dust*/
	  r          = tau/(lya_sigma*N + dust_sigma*All.NumberDensityDust); /*the distance is in cm*/
	}

	r_travel[i_photon] = r;

	if((rand_interaction > HIInteractionProb) && All.UseDust && (rand_absorption < All.DustAbsorptionProb)){
	  r_travel[i_photon] = 0.0;
	  status[i_photon] = ABSORBED;
	  x_out[i_photon] = x_in[i_photon];
	}
      }
    }
    return 0;
}

void OpenAsciiFile(char *fname)
{
    FILE *f;
    if(!(f = fopen(fname,"w")))
    {
	fprintf(stderr,"DumpPhotonList: Problem opening file %s\n",fname);
	exit(0);

    }
    fprintf(f, "# %d %e %e %e %e %e\n", 50000, All.Tau, All.Temperature, All.InputFrequency, All.TauDust, All.DustAbsorptionProb);
    fclose(f);
}


void AppendAsciiFile(char *fname, double PosX, double PosY, double PosZ, double DirX, double DirY, double DirZ, double x_out, int status, int n_scatter)
{
    FILE *f;

    if(!(f = fopen(fname,"a")))
    {
	fprintf(stderr,"DumpPhotonList: Problem opening file %s\n",fname);
	exit(0);
    }

    fprintf(f,"%e %e %e %e %e %e %e %d %d\n", 
	    PosX, PosY, PosZ,
	    DirX, DirY, DirZ, 
	    x_out, status, n_scatter);
    fclose(f);
}


extern "C" int PropagatePackage(double *PosX, double *PosY, double *PosZ, 
		     double *DirX, double *DirY, double *DirZ, int *n_scatt, 
		     double *x_in, int *status, int n_points){
    int i;
    double Pos[3];
    double Dir[3];
    double r_travel, x, x_comoving;
    double a, nu_doppler, n_HI, v_thermal, BulkVel[3], temperature;
    int n_iter;
    int stat;
    float last_x;
    int  n_active;
    int n_global_scatt;
    float *x_aux_in;
    float *x_aux_out;
    double *r_travel_aux;
    float *v_parallel;
    float *v_perp_1;
    float *v_perp_2;
    float *v_perp_1_d;
    float *v_perp_2_d;
    float *v_parallel_d;
    float *x_aux_d;
    curandState *d_rngStates = 0;
    cudaError_t cudaResult = cudaSuccess;
    int nBlocks, blockSize;
    n_iter=0;
    n_global_scatt=0;
    char FileName[MAX_FILENAME_SIZE];
    int n_total;

    if(All.OutputFinalList){
	sprintf(FileName, "%s/%s_out.ascii", All.OutputDir, All.OutputFile);
	OpenAsciiFile(FileName);
    }

    /* auxiliary variables */
    if(!(x_aux_in = (float *)malloc(sizeof(float) * n_points))){
      fprintf(stderr, "Problem with x_aux allocation\n");
      exit(1);
    }

    if(!(x_aux_out = (float *)malloc(sizeof(float) * n_points))){
      fprintf(stderr, "Problem with x_aux allocation\n");
      exit(1);
    }

    if(!(r_travel_aux = (double *)malloc(sizeof(double) * n_points))){
      fprintf(stderr, "Problem with r_travel_aux allocation\n");
      exit(1);
    }

    if(!(v_parallel = (float *)malloc(sizeof(float) * n_points))){
      fprintf(stderr, "Problem with r_travel_aux allocation\n");
      exit(1);
    }

    if(!(v_perp_1 = (float *)malloc(sizeof(float) * n_points))){
      fprintf(stderr, "Problem with r_travel_aux allocation\n");
      exit(1);
    }

    if(!(v_perp_2 = (float *)malloc(sizeof(float) * n_points))){
      fprintf(stderr, "Problem with r_travel_aux allocation\n");
      exit(1);
    }

    /*Device Malloc*/
    cudaMalloc((void **) &v_parallel_d, n_points * sizeof(float));
    cudaMalloc((void **) &x_aux_d, n_points * sizeof(float));


    
    blockSize = 128; // This is the number of threads inside a block
    nBlocks = (n_points)/blockSize + (n_points%blockSize == 0?0:1); // This is the number of blocks
    fprintf(stdout, "nBlocks %d\n", nBlocks);
    
    /*RNG*/
    cudaResult = cudaMalloc((void **)&d_rngStates, blockSize * nBlocks * sizeof(curandState));
    initRNG<<<nBlocks, blockSize>>>(d_rngStates, All.RandomSeed);
   




    n_total = 0;
    while(n_total<5000){
    /*Get all the atom velocities -  This is the part to be updated with a kernel*/
    for(i=0;i<n_points;i++){
      /*If the photon is not active anymore, I re-initialize its values*/
      if(status[i]!=ACTIVE){
#ifdef DEBUG
	fprintf(stdout, "Total photons: %d\n", n_total);
	fflush(stdout);
#endif
	n_total++;
	AppendAsciiFile(FileName, PosX[i], PosY[i], PosZ[i], DirX[i], DirY[i], DirZ[i], x_in[i], status[i], n_scatt[i]);
	x_in[i] = 0.0;
	n_scatt[i] = 0;
	PhotonInitialize(&(PosX[i]), &(PosY[i]), &(PosZ[i]), &(DirX[i]), &(DirY[i]), &(DirZ[i]));
	status[i] = ACTIVE;
      }

      Pos[0] = PosX[i];
      Pos[1] = PosY[i];
      Pos[2] = PosZ[i];
      stat = status[i];
      x_aux_in[i] = x_in[i];    
      x_aux_out[i] = x_in[i];    
      
      /* get the temperature at this point*/
      PropagateGetTemperature(&temperature, Pos);
      
      /*Get the thermal velocity and doppler broadening*/
      nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
      a = Lya_nu_line_width_CGS/(2.0*nu_doppler);		
    }
      
    /****THIS IS THE GPU PART*****/
    /*get first the parallel velocity*/
    cudaMemcpy(x_aux_d, x_aux_in, sizeof(float) * n_points, cudaMemcpyHostToDevice);
    get_lyman_parallel_vel<<<nBlocks, blockSize>>>(v_parallel_d, x_aux_d, a, d_rngStates, n_points);      
    cudaMemcpy(v_parallel, v_parallel_d, sizeof(float) * n_points, cudaMemcpyDeviceToHost);
    
    /*get the perpendicular velocity*/
    get_lyman_perp_vel(v_perp_1, v_perp_2, n_points);	        
    /***************************/
    
    for(i=0;i<n_points;i++){
      /*Make the initialization*/
      Pos[0] = PosX[i];
      Pos[1] = PosY[i];
      Pos[2] = PosZ[i];
      Dir[0] = DirX[i];
      Dir[1] = DirY[i];
      Dir[2] = DirZ[i];
      stat = status[i];
      x_aux_in[i] = x_in[i];    
      x_aux_out[i] = x_in[i];    
      
      /*If the photon is stil inside and active, update its direction of propagation and frequency 
	to  gas rest frame */
      if(PropagateIsInside(Pos[0],Pos[1],Pos[2])){
	/* get the temperature at this point*/
	PropagateGetTemperature(&temperature, Pos);
	
	/* get the number density at this point*/
	PropagateGetNumberDensity(&n_HI, Pos);
	
	/*get the bulk velocity of the fluid at this point*/
	PropagateGetBulkVel(BulkVel, Pos);
	
	/*Get the thermal velocity and doppler broadening*/
	nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
	a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
	v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
	
	/*change the value of the frequency to one comoving with the fluid*/
	PropagateLorentzFreqChange(&(x_aux_in[i]), Dir, BulkVel, v_thermal, -1); 
	
	/*change the direction of the photon to the fluid frame*/
	PropagateLorentzDirChange(&(Dir[0]), BulkVel, -1);
      }
      DirX[i] = Dir[0];
      DirY[i] = Dir[1];
      DirZ[i] = Dir[2];
    }
    
    PropagateGetNumberDensity(&n_HI, Pos);
    PropagateGetTemperature(&temperature, Pos);
    nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
    a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
    v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
    
    /*Change the frequency and the Propagation direction, find the displacement*/	
    PropagateStep(v_parallel, v_perp_1, v_perp_2, x_aux_in, x_aux_out, DirX, DirY, DirZ, r_travel_aux, status, a, n_HI, n_points);	    	
    
    for(i=0;i<n_points;i++){
      Pos[0] = PosX[i];
      Pos[1] = PosY[i];
      Pos[2] = PosZ[i];	
      Dir[0] = DirX[i];
      Dir[1] = DirY[i];
      Dir[2] = DirZ[i];
      stat = status[i];	

      if(PropagateIsInside(Pos[0],Pos[1],Pos[2])&&(stat==ACTIVE)){
	/* get the temperature at this point*/
	PropagateGetTemperature(&temperature, Pos);
	
	/* get the number density at this point*/
	PropagateGetNumberDensity(&n_HI, Pos);
	  
	  /*get the bulk velocity of the fluid at this point*/
	  PropagateGetBulkVel(BulkVel, Pos);

	  /*Get the thermal velocity and doppler broadening*/
	  nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
	  a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
	  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
	  
	  /*Change the new direction to the lab frame value*/
	  PropagateLorentzDirChange(Dir, BulkVel, 1);
	  
	  /*Change the frequency comoving to the lab frame value*/
	  PropagateLorentzFreqChange(&(x_aux_out[i]), Dir, BulkVel, v_thermal, 1); 
	  
	  /*Update the position*/
	  Pos[0] += r_travel_aux[i] * Dir[0];	    
	  Pos[1] += r_travel_aux[i] * Dir[1];	    
	  Pos[2] += r_travel_aux[i] * Dir[2];    
	  n_scatt[i]++;


	  /*update the final status of the photon, just to know if it was absorbed, or what*/
	  if(!(PropagateIsInside(Pos[0],Pos[1],Pos[2]))){
	    status[i] = OUT_OF_BOX;
	  }

	  PosX[i] = Pos[0];
	  PosY[i] = Pos[1];
	  PosZ[i] = Pos[2];
	  DirX[i] = Dir[0];
	  DirY[i] = Dir[1];
	  DirZ[i] = Dir[2];	  
	  x_in[i] = x_aux_out[i];
      }
    }
    }
    

    free(x_aux_in);
    free(x_aux_out);
    free(r_travel_aux);
    free(v_parallel);
    free(v_perp_1);
    free(v_perp_2);
    return 0;    
}
