#ifndef PROP_KERNEL_H
#define PROP_KERNEL_H

extern "C" 
int PropagatePackage(double *PosX, double *PosY, double *PosZ, 
				double *DirX, double *DirY, double *DirZ, int *n_scatt, 
				double *x_in, int *status, int n_points);

#endif
