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

