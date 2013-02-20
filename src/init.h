#ifndef INIT_H
#define INIT_H
float * InitCreationFrequency(int n_points);
float * InitCreationPosition(int n_points);
float * InitCreationDirection(int n_points);

void InitFrequency(float *x, int n_points);
void InitPosition(float *p, int n_points);
void InitDirection(float *k, int n_points);
void InitCheckDirection(float *k, int n_points);
#endif
