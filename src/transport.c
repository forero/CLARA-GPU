#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "transport.h"
#include "scatter.h"

void TransportPhotons(float *x, int n_photons){
  int pack_size, last_pack_size;  
  int n_packs;
  int i;

  pack_size = 400;
  
  n_packs = n_photons/pack_size;
  last_pack_size = n_photons%pack_size;
  
  fprintf(stdout, "Photons to transport %d\n", n_photons);
  fprintf(stdout, "%d packs of size %d\n", n_packs, pack_size);
  fprintf(stdout, "one last pack of size %d\n", last_pack_size);
  
  for(i=0;i<n_packs;i++){
    scatter_x(x, i*pack_size, (i+1)*pack_size);
  }
}
