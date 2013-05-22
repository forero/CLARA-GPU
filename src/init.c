#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "init.h"

int * InitCreationInteger(int n_points){
  int *i;
  if(!(i=malloc(n_points * sizeof(int)))){
    fprintf(stderr, "Problem in the allocation of integer array\n");
    exit(1);
  }
  return i;
}


FLOAT * InitCreationFrequency(int n_points){
  FLOAT *x;
  if(!(x=malloc(n_points * sizeof(FLOAT)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
    return x;
}


FLOAT * InitCreationPosition(int n_points){
  FLOAT *p;
  if(!(p=malloc(3 * n_points * sizeof(FLOAT)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
  return p;
}

FLOAT * InitCreationDirection(int n_points){
  FLOAT *k;
  if(!(k=malloc(3 * n_points * sizeof(FLOAT)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
  return k;
}


void InitNScatt(int *n, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    n[i] = 0;
  }
}

void InitStatus(int *n, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    n[i] = ACTIVE;
  }
}

void InitFrequency(FLOAT *x, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    x[i] = 0.0;
  }
}

void InitPosition(FLOAT *p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[3*i + 0] = 0.0;
    p[3*i + 1] = 0.0;
    p[3*i + 2] = 0.0;
  }
}


void InitDirection(FLOAT *k, int n_points){
  int i;
  FLOAT theta, phi;
  /*random directions over the sphere*/
  for(i=0;i<n_points;i++){
    theta = acos((drand48()-0.5)*2.0);
    phi = 2.0*PI*drand48();    
    k[3*i + 0] = sin(theta)*cos(phi);
    k[3*i + 1] = sin(theta)*sin(phi);
    k[3*i + 2] = cos(theta);
  }  
}

void InitCheckDirection(FLOAT *k, int n_points){
  int i;
  FLOAT norm2;
  for(i=0;i<n_points;i++){
      norm2 = 
	(k[3*i + 0]*k[3*i + 0])  +
	(k[3*i + 1] * k[3*i + 1]) +
	(k[3*i + 2] * k[3*i + 2]);
      if(fabs(norm2-1.0)>EPSILON_DOUBLE){	
	fprintf(stderr, "Directions are not normalized to unity. EXIT.\n");
	fprintf(stderr, "photon id %d. %f %f %f\n", i, k[3*i + 0],k[3*i + 1],k[3*i + 2]);
	exit(1);
      }
  }  

}
