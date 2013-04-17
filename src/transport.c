#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include "struct.h"
#include "scatter.h"


void TransportPhotons(FLOAT *x, FLOAT *p, FLOAT *k, int n_photons){
  int pack_size, last_pack_size;  
  int n_packs;
  int i;

  pack_size = 512;
  
  n_packs = n_photons/pack_size;
  last_pack_size = n_photons%pack_size;
  
  fprintf(stdout, "Photons to transport %d\n", n_photons);
  fprintf(stdout, "%d packs of size %d\n", n_packs, pack_size);
  fprintf(stdout, "one last pack of size %d\n", last_pack_size);
  
  for(i=0;i<n_packs;i++){
    scatter_bunch(x, p, k, i*pack_size, (i+1)*pack_size);
  }  
  scatter_bunch(x, p, k, i*pack_size, i*pack_size + last_pack_size);  
}
