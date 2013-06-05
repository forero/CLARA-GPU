#include <stdio.h>
#include <math.h>

__device__ int PropagateIsInside(FLOAT pos_x, FLOAT pos_y, FLOAT pos_z, setup *S)
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

__device__ void PropagateGetTemperature(FLOAT *Temp, FLOAT *Pos, setup *S){
    *Temp = S->Temperature;
    return;
}

__device__ void PropagateGetNumberDensity(FLOAT *n_HI, FLOAT *Pos, setup *S){
    *n_HI = S->NumberDensityHI;
}

__device__ void PropagateLorentzDirChange(FLOAT *Dir, FLOAT *BulkVel, int sign){
    int i;
    for(i=0;i<3;i++){
	Dir[i] = Dir[i]*(1.0 + sign*(BulkVel[i]/C_LIGHT)); 
    }    
    normalize(Dir);
}

__device__ void PropagateLorentzFreqChange(FLOAT *x, FLOAT *Dir, 
			    FLOAT *BulkVel, FLOAT ThermalVel, int sign){
    FLOAT lorentz_factor;
    lorentz_factor = (BulkVel[0]*Dir[0] + BulkVel[1]*Dir[1] + BulkVel[2]*Dir[2])/ThermalVel;
    
    *x  = *x + sign*lorentz_factor;
}

__device__ void PropagateGetBulkVel(FLOAT *BulkVel, FLOAT *Pos, setup *S)
/*Returns bulk velocity in cm/s*/
{
    int i;

    for(i=0;i<3;i++){
	BulkVel[i] = 0.0;
    }

    if(S->ExpandingSphere){
      for(i=0;i<3;i++){
	BulkVel[i] = (Pos[i]/S->Tau) * S->Vmax * 1.0e5;
      }	
    }

    if(S->RotatingSphere){
      BulkVel[0]=(-Pos[1]/S->Tau) * S->Vmax * 1.0e5;
      BulkVel[1]=(Pos[0]/S->Tau) * S->Vmax * 1.0e5;
      BulkVel[2]=0.0;
    }	
}


__device__ FLOAT LyaTau(curandState *state)
/*generates a value of optical depth*/
{
    FLOAT tau;
    FLOAT tmp;
    tmp = curand_uniform(state);
    tau = -log(tmp);
    return tau;
}


__device__ FLOAT LyaH(FLOAT x, FLOAT a)
{
    FLOAT z, P, q, H;
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


__device__ FLOAT LyaCrossSection(FLOAT x, FLOAT a){
    FLOAT sigma;
    FLOAT nu_doppler;


    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    sigma = 0.4162*sqrt(PI)*(CHARGEELECTRON*CHARGEELECTRON)/(ELECTRONMASS*C_LIGHT*nu_doppler);
    sigma = sigma*LyaH(x,a);
        
    return sigma;
}


__device__ int PropagateStep(FLOAT *x, FLOAT *Dir, FLOAT *tau_travel, FLOAT *a, 
			     curandState *state, setup *S)
/*
  This fuction calculates the new frequency, direction and traveled
 distance by a single photon
 
 Important: at this point the frequency must be the comoving one (not lab frame!). If the photon is absorbed by dust the program returns -1
*/
{
  FLOAT nu_doppler;
  FLOAT v_thermal;
  FLOAT temperature;
  FLOAT g_recoil;
  FLOAT travel_tau;
  FLOAT lya_sigma_0, lya_sigma, lya_sigma_ratio;

  /*variable to define the interaction with an atom or a 
   dust grain*/  
  FLOAT HIInteractionProb;
  FLOAT rand_absorption;
  FLOAT rand_interaction;

  FLOAT x_in, x_out, a_in;
  int status;
  int program_status;
  int i;
  /*atom properties*/
  FLOAT u_atom[3];
  FLOAT k_out_photon[3];
  FLOAT k_in_photon[3];
  

  /*variable initialization*/
  status = ACTIVE;
  x_in = *x;
  a_in = *a;
  for(i=0;i<3;i++){
    k_in_photon[i] = Dir[i];
  }

  /*recompute important physical quantities in cgs units*/
  nu_doppler = Lya_nu_line_width_CGS/(2.0*a_in);
  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
  temperature = (nu_doppler/CONSTANT_NU_DOPPLER)*
    (nu_doppler/CONSTANT_NU_DOPPLER)*10000.0;
  g_recoil = PLANCK*nu_doppler/(2.0*BOLTZMANN*temperature);
  
  /*gas status at the interaction point*/
  lya_sigma = LyaCrossSection(x_in, a_in);
  lya_sigma_0 = LyaCrossSection(0, a_in);
  lya_sigma_ratio = lya_sigma / lya_sigma_0;
  travel_tau = LyaTau(state);

  /*get a random number to decide if the thing interact with HI*/
  rand_interaction  = curand_uniform(state);
  HIInteractionProb = 1.0;
  rand_absorption = 0.0;

  
  HIInteractionProb = (S->Tau*lya_sigma_ratio)/((S->Tau*lya_sigma_ratio) + S->TauDust);
  rand_absorption = curand_uniform(state);


  /*If present, dust changes the traveled tau*/
  travel_tau = travel_tau * HIInteractionProb;

  
  if(rand_interaction <= HIInteractionProb){
    RND_lyman_atom(&(u_atom[0]), &(k_in_photon[0]), &(k_out_photon[0]), x_in, a_in , state, &program_status);
    
    /*find the new frequency (in the observer system)*/
    x_out = x_in - point_product(u_atom, k_in_photon) +
      point_product(u_atom, k_out_photon) +	  
      g_recoil * (point_product(k_in_photon, k_out_photon) - 1.0);

  }else{
    x_out = x_in; 
  }

  if((rand_interaction > HIInteractionProb) && (S->TauDust>0.0) && (rand_absorption < S->DustAbsorptionProb)){
    travel_tau = 0.0;
    status = ABSORBED;
  }
  
  
  /*Update the relevant physical variables*/
  *x = x_out;
  *tau_travel = travel_tau;
  for(i=0;i<3;i++){
    Dir[i] = k_out_photon[i];
  }  
  
  return status;
}
