#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "init.h"
#include "scatter.h"
#include "io.h"
#include "transport.h"
#include "parse.h"
#include "struct.h"

int main(int argc, char **argv){
  float *x;
  float *p;
  float *k;
  int n_points;
    
  /*parse the parameter setup*/
  ReadParameters(argv[1]);

  /*Create the photons*/
  x = InitCreationFrequency(All.NPackages);
  p = InitCreationPosition(All.NPackages);
  k = InitCreationDirection(All.NPackages);

  /*Initialize them*/
  InitFrequency(x, All.NPackages);
  InitPosition(p, All.NPackages);
  InitDirection(k, All.NPackages);

  /*Write files to disk*/
  DumpPhotonList(x, p, k, All.NPackages, "test.in");

  /*Sanity Checks Before RT*/
  InitCheckDirection(k, All.NPackages);
  
  /*Make the radiative transfer*/
  TransportPhotons(x, p, k, All.NPackages);

  /*Sanity Checks after RT*/
  InitCheckDirection(k, All.NPackages);

  /*Write files to disk*/
  DumpPhotonList(x, p, k, All.NPackages, "test.out");
  
  return 0;
}

