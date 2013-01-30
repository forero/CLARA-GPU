#include <stdio.h>
#include <stdlib.h>
#include "init.h"

float * InitCreation(int n_points){
  float *x;
  if(!(x=malloc(n_points * sizeof(float)))){
    fprintf(stderr, "Problem in the allocation of property x\n");
    exit(1);
  }
    return x;
}

void InitFrequency(float *x, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    x[i] = 0.0;
  }
}
