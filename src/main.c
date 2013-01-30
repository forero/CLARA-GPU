#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "init.h"
#include "scatter.h"

int main(int argc, char **argv){
  float *x;
  int n_points;

  n_points = atoi(argv[1]);

  x = InitCreation(n_points);

  InitFrequency(x, n_points);
  
  TransportPhotons(x, n_points);
  
  return 0;
}

