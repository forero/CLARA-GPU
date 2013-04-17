#include "struct.h"
#ifndef INIT_H
#define INIT_H
FLOAT * InitCreationFrequency(int n_points);
FLOAT * InitCreationPosition(int n_points);
FLOAT * InitCreationDirection(int n_points);

void InitFrequency(FLOAT *x, int n_points);
void InitPosition(FLOAT *p, int n_points);
void InitDirection(FLOAT *k, int n_points);
void InitCheckDirection(FLOAT *k, int n_points);
#endif
