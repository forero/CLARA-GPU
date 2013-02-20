#include <stdio.h>
#include <stdlib.h>

void DumpPhotonList(float *x, float *p, float *k, int n_points, char *filename){
  FILE *out;
  int i;

  if(!(out=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  
  for(i=0;i<n_points;i++){
    fprintf(out, "%f %f %f %f %f %f %f\n", x[i], 
	    p[3*i + 0], p[3*i + 1], p[3*i +2], 
	    k[3*i + 0], k[3*i + 1], k[3*i +2]);	    
  }
  fclose(out);
}
