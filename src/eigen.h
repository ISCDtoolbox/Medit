#ifndef _EIGEN_H
#define _EIGEN_H

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#define EV_EPS    1.e-6
#define EV_EPS1   1.e-13
#define EV_EPSD   1.e-200

enum {None=0, SymMat=1};


/* check if numbers are equal */ 
#define egal(x,y)   ( \
  ( ((x) == 0.0) ? (fabs(y) < EV_EPS) : \
  ( ((y) == 0.0) ? (fabs(x) < EV_EPS) : \
  (fabs((x)-(y)) / (fabs(x) + fabs(y)) < EV_EPS) ) ) )


/* prototypes */
int eigen_3d(int symmat,double *mat,double lambda[3],double v[3][3]);
int eigen_2d(double *m,double *l,double vp[2][2]);

#endif