#include <stdio.h>
#include <stdlib.h>

void DumpPhotonList(float *list, int n_points, char *filename){
  FILE *out;
  int i;

  if(!(out=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  
  for(i=0;i<n_points;i++){
    fprintf(out, "%f\n", list[i]);
  }
  fclose(out);
}
