#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init.h"
#include "struct.h"

float * InitCreationFrequency(int n_points){
  float *x;
  if(!(x=malloc(n_points * sizeof(float)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
    return x;
}


float * InitCreationPosition(int n_points){
  float *p;
  if(!(p=malloc(3 * n_points * sizeof(float)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
  return p;
}

float * InitCreationDirection(int n_points){
  float *k;
  if(!(k=malloc(3 * n_points * sizeof(float)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
  return k;
}


void InitFrequency(float *x, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    x[i] = 0.0;
  }
}

void InitPosition(float *p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[3*i + 0] = 0.0;
    p[3*i + 1] = 0.0;
    p[3*i + 2] = 0.0;
  }
}


void InitDirection(float *k, int n_points){
  int i;
  float theta, phi;
  /*random directions over the sphere*/
  for(i=0;i<n_points;i++){
    theta = acos((drand48()-0.5)*2.0);
    phi = 2.0*PI*drand48();    
    k[3*i + 0] = sin(theta)*cos(phi);
    k[3*i + 1] = sin(theta)*sin(phi);
    k[3*i + 2] = cos(theta);
  }  
}

void InitCheckDirection(float *k, int n_points){
  int i;
  float norm2;
  for(i=0;i<n_points;i++){
      norm2 = (k[3*i + 0]*k[3*i + 0])  +
	(k[3*i + 1] * k[3*i + 1]) +
	(k[3*i + 2] * k[3*i + 2]);
      if(fabs(norm2-1.0)>EPSILON_DOUBLE){
	fprintf(stderr, "Directions are not normalized to unity. EXIT.\n");
	exit(1);
      }
  }  

}
