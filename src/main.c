#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "init.h"
#include "scatter.h"
#include "io.h"

int main(int argc, char **argv){
  float *x;
  float *p;
  float *k;
  int n_points;

  /*get the parameter setup*/
  n_points = atoi(argv[1]);

  /*Create the photons*/
  x = InitCreationFrequency(n_points);
  p = InitCreationPosition(n_points);
  k = InitCreationDirection(n_points);

  /*Initialize them*/
  InitFrequency(x, n_points);
  InitPosition(p, n_points);
  InitDirection(k, n_points);

  /*Sanity Checks Before RT*/
  InitCheckDirection(k, n_points);
  
  /*Make the radiative transfer*/
  TransportPhotons(x, p, k, n_points);

  /*Sanity Checks after RT*/
  InitCheckDirection(k, n_points);

  /*Write files to disk*/
  DumpPhotonList(x, p, k, n_points, "test.out");
  
  return 0;
}

