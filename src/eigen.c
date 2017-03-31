#include "eigen.h"


static double Id[3][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
#define EV_MAXIT    50


/* find root(s) of polynomial:  P(x)= x^3+bx^2+cx+d 
   return 1: 3 roots,  2: 2 roots,  3: 1 root */
static int newton3(double p[4],double x[3]) {
  double     b,c,d,da,db,dc,epsd;
  double     delta,fx,dfx,dxx;
  double     fdx0,fdx1,dx0,dx1,x1,x2;
  int        it,n;

  /* coeffs polynomial, a=1 */
  b = p[2];
  c = p[1];
  d = p[0];
  n = 1;

  /* 1st derivative of f */
  da = 3.0;
  db = 2.0*b;

  /* solve 2nd order eqn */
  delta = db*db - 4.0*da*c;
  epsd  = db*db*EV_EPSD;

  /* inflexion (f'(x)=0, x=-b/2a) */
  x1 = -db / 6.0;

  n = 1;
  if ( delta > epsd ) {
    delta = sqrt(delta);
    dx0   = (-db + delta) / 6.0;
    dx1   = (-db - delta) / 6.0;
    /* Horner */
    fdx0 = d + dx0*(c+dx0*(b+dx0));
    fdx1 = d + dx1*(c+dx1*(b+dx1));

    if ( fabs(fdx0) < EV_EPSD ) {
      /* dx0: double root, compute single root */
      n = 2;
      x[0] = dx0;
      x[1] = dx0;
      x[2] = -b - 2.0*dx0;
      /* check if P(x) = 0 */
      fx = d + x[2]*(c+x[2]*(b+x[2]));
      if ( fabs(fx) > EV_EPSD )
        return(0);
      else
        return(n);
    }
    else if ( fabs(fdx1) < EV_EPSD ) {
      /* dx1: double root, compute single root */
      n = 2;
      x[0] = dx1;
      x[1] = dx1;
      x[2] = -b - 2.0*dx1;
      /* check if P(x) = 0 */
      fx = d + x[2]*(c+x[2]*(b+x[2]));
      if ( fabs(fx) > EV_EPSD )
        return(0);
      else
        return(n);
    }
  }
  else if ( fabs(delta) < epsd ) {
    /* triple root */
    n = 3;
    x[0] = x1;
    x[1] = x1;
    x[2] = x1;
    /* check if P(x) = 0 */
    fx = d + x[0]*(c+x[0]*(b+x[0]));
    if ( fabs(fx) > EV_EPSD )
      return(0);
    else
      return(n);
  }
  else {
    return(0);
  }

  /* Newton method: find one root (middle)
     starting point: P"(x)=0 */
  x1  = -b / 3.0;
  dfx =  c + b*x1;
  fx  = d + x1*(c -2.0*x1*x1);
  it  = 0;
  do {
    x2 = x1 - fx / dfx;
    fx = d + x2*(c+x2*(b+x2));
    if ( fabs(fx) < EV_EPSD ) {
      x[0] = x2;
      break;
    }
    dfx = c + x2*(db + da*x2);

    /* check for break-off condition */
    dxx = fabs((x2-x1) / x2);
    if ( dxx < 1.0e-10 ) {
      x[0] = x2;
      if ( fabs(fx) > EV_EPSD )  return(0);
      break;
    }
    else
      x1 = x2;
  }
  while ( ++it < EV_MAXIT );

  if ( it == EV_MAXIT ) {
    x[0] = x1;
    fx   = d + x1*(c+(x1*(b+x1)));
    if ( fabs(fx) > EV_EPSD )  return(0);
  }

  /* solve 2nd order equation
     P(x) = (x-sol(1))* (x^2+bb*x+cc)  */
  db    = b + x[0];
  dc    = c + x[0]*db;
  delta = db*db - 4.0*dc;
  
  if ( delta <= 0.0 )  return(0);

  delta = sqrt(delta);
  x[1] = 0.5 * (-db+delta);
  x[2] = 0.5 * (-db-delta);

  return(n);
}


/* find eigenvalues and vectors of a 3x3 symmetric definite
 * positive matrix. return order of eigenvalues (1,2,3) or 0 */
int eigen_3d(int type,double *m,double *l,double v[3][3]) {
  double    a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double    aa,bb,cc,dd,ee,ii,vx1[3],vx2[3],vx3[3],dd1,dd2,dd3;
  double    maxd,maxm,valm,p[4],w1[3],w2[3],w3[3];
  int       k,n;

  /* default */
  memcpy(v,Id,9*sizeof(double));
  if ( type == SymMat ) {
    l[0] = m[0];
    l[1] = m[3];
    l[2] = m[5];

    maxm = fabs(m[0]);
    for (k=1; k<6; k++) {
      valm = fabs(m[k]);
      if ( valm > maxm )  maxm = valm;
    }
    if ( maxm > EV_EPSD ) {
      /* normalize matrix */
      a11 = m[0] / maxm;
      a12 = m[1] / maxm;
      a13 = m[2] / maxm;
      a22 = m[3] / maxm;
      a23 = m[4] / maxm;
      a33 = m[5] / maxm;
    }
    /* diagonal matrix */
    maxd = fabs(a12);
    valm = fabs(a13);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a23);
    if ( valm > maxd )  maxd = valm;
    if ( maxd < EV_EPSD )  return(1);

    a21  = a12;
    a31  = a13;
    a32  = a23;

    /* build characteristic polynomial
       P(X) = X^3 - trace X^2 + (somme des mineurs)X - det = 0 */
    aa = a11*a22;
    bb = a23*a32;
    cc = a12*a21;
    dd = a13*a31;
    p[0] =  a11*bb + a33*(cc-aa) + a22*dd -2.0*a12*a13*a23;
    p[1] =  a11*(a22 + a33) + a22*a33 - bb - cc - dd;
    p[2] = -a11 - a22 - a33;
    p[3] =  1.0;
  }
  else {
    l[0] = m[0];
    l[1] = m[4];
    l[2] = m[8];

    maxm = fabs(m[0]);
    for (k=1; k<9; k++) {
      valm = fabs(m[k]);
      if ( valm > maxm )  maxm = valm;
    }
    if ( maxm > EV_EPSD ) {
      /* normalize matrix */
      a11 = m[0] / maxm;
      a12 = m[1] / maxm;
      a13 = m[2] / maxm;
      a21 = m[3] / maxm;
      a22 = m[4] / maxm;
      a23 = m[5] / maxm;
      a31 = m[6] / maxm;
      a32 = m[7] / maxm;
      a33 = m[8] / maxm;
    }

    /* diagonal matrix */
    maxd = fabs(a12);
    valm = fabs(a13);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a23);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a21);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a31);
    if ( valm > maxd )  maxd = valm;
     valm = fabs(a32);
    if ( valm > maxd )  maxd = valm;
    if ( maxd < EV_EPSD )  return(1);

    /* build characteristic polynomial
       P(X) = X^3 - trace X^2 + (somme des mineurs)X - det = 0 */
    aa = a22*a33 - a23*a32;
    bb = a23*a31 - a21*a33;
    cc = a21*a32 - a31*a22;
    ee = a11*a33 - a13*a31;
    ii = a11*a22 - a12*a21;
    
    p[0] =  -a11*aa - a12*bb - a13*cc;
    p[1] =  aa + ee + ii;
    p[2] = -a11 - a22 - a33;
    p[3] =  1.0;
  }

  /* solve polynomial (find roots using newton) */
  n = newton3(p,l);
  if ( n <= 0 )  return(0);

  /* compute eigenvectors:
     an eigenvalue belong to orthogonal of Im(A-lambda*Id) */
  v[0][0] = 1.0; v[0][1] = v[0][2] = 0.0;
  v[1][1] = 1.0; v[1][0] = v[1][2] = 0.0;
  v[2][2] = 1.0; v[2][0] = v[2][1] = 0.0;

  w1[1] = a12;  w1[2] = a13;
  w2[0] = a21;  w2[2] = a23;
  w3[0] = a31;  w3[1] = a32;

  if ( n == 1 ) {
    /* vk = crsprd(wi,wj) */
    for (k=0; k<3; k++) {
      w1[0] = a11 - l[k];
      w2[1] = a22 - l[k];
      w3[2] = a33 - l[k];

      /* cross product vectors in (Im(A-lambda(i) Id) ortho */
      vx1[0] = w1[1]*w3[2] - w1[2]*w3[1];
      vx1[1] = w1[2]*w3[0] - w1[0]*w3[2];
      vx1[2] = w1[0]*w3[1] - w1[1]*w3[0];
      dd1    = vx1[0]*vx1[0] + vx1[1]*vx1[1] + vx1[2]*vx1[2];

      vx2[0] = w1[1]*w2[2] - w1[2]*w2[1];
      vx2[1] = w1[2]*w2[0] - w1[0]*w2[2];
      vx2[2] = w1[0]*w2[1] - w1[1]*w2[0];
      dd2    = vx2[0]*vx2[0] + vx2[1]*vx2[1] + vx2[2]*vx2[2];

      vx3[0] = w2[1]*w3[2] - w2[2]*w3[1];
      vx3[1] = w2[2]*w3[0] - w2[0]*w3[2];
      vx3[2] = w2[0]*w3[1] - w2[1]*w3[0];
      dd3    = vx3[0]*vx3[0] + vx3[1]*vx3[1] + vx3[2]*vx3[2];

      /* find vector of max norm */
      if ( dd1 > dd2 ) {
        if ( dd1 > dd3 ) {
          dd1 = 1.0 / sqrt(dd1);
          v[k][0] = vx1[0] * dd1;
          v[k][1] = vx1[1] * dd1;
          v[k][2] = vx1[2] * dd1;
        }
	      else {
          dd3 = 1.0 / sqrt(dd3);
          v[k][0] = vx3[0] * dd3;
          v[k][1] = vx3[1] * dd3;
          v[k][2] = vx3[2] * dd3;
	      }
      }
      else {
        if ( dd2 > dd3 ) { 
          dd2 = 1.0 / sqrt(dd2);
          v[k][0] = vx2[0] * dd2;
          v[k][1] = vx2[1] * dd2;
          v[k][2] = vx2[2] * dd2;
        }
        else {
          dd3 = 1.0 / sqrt(dd3);
          v[k][0] = vx3[0] * dd3;
          v[k][1] = vx3[1] * dd3;
          v[k][2] = vx3[2] * dd3;
        }  
      }
    }
  }

  /* (vp1,vp2) double,  vp3 simple root */
  else if ( n == 2 ) {
    w1[0] = a11 - l[2];
    w2[1] = a22 - l[2];
    w3[2] = a33 - l[2];

    /* cross product */
    vx1[0] = w1[1]*w3[2] - w1[2]*w3[1];
    vx1[1] = w1[2]*w3[0] - w1[0]*w3[2];
    vx1[2] = w1[0]*w3[1] - w1[1]*w3[0];
    dd1 = vx1[0]*vx1[0] + vx1[1]*vx1[1] + vx1[2]*vx1[2];
 
    vx2[0] = w1[1]*w2[2] - w1[2]*w2[1];
    vx2[1] = w1[2]*w2[0] - w1[0]*w2[2];
    vx2[2] = w1[0]*w2[1] - w1[1]*w2[0];
    dd2 = vx2[0]*vx2[0] + vx2[1]*vx2[1] + vx2[2]*vx2[2];

    vx3[0] = w2[1]*w3[2] - w2[2]*w3[1];
    vx3[1] = w2[2]*w3[0] - w2[0]*w3[2];
    vx3[2] = w2[0]*w3[1] - w2[1]*w3[0];
    dd3 = vx3[0]*vx3[0] + vx3[1]*vx3[1] + vx3[2]*vx3[2];

    /* find vector of max norm */
    if ( dd1 > dd2 ) {
      if ( dd1 > dd3 ) {
        dd1 = 1.0 / sqrt(dd1);
        v[2][0] = vx1[0] * dd1;
        v[2][1] = vx1[1] * dd1;
        v[2][2] = vx1[2] * dd1;
      }
      else {
        dd3 = 1.0 / sqrt(dd3);
        v[2][0] = vx3[0] * dd3;
        v[2][1] = vx3[1] * dd3;
        v[2][2] = vx3[2] * dd3;
      }
    }
    else {
      if ( dd2 > dd3 ) {
        dd2 = 1.0 / sqrt(dd2);
        v[2][0] = vx2[0] * dd2;
        v[2][1] = vx2[1] * dd2;
        v[2][2] = vx2[2] * dd2;
      }
      else {
        dd3 = 1.0 / sqrt(dd3);
        v[2][0] = vx3[0] * dd3;
        v[2][1] = vx3[1] * dd3;
        v[2][2] = vx3[2] * dd3;
      }
    }

    /* compute v1 and v2 in Im(A-vp3*Id) */
    dd1 = w1[0]*w1[0] + w1[1]*w1[1] + w1[2]*w1[2];
    dd2 = w2[0]*w2[0] + w2[1]*w2[1] + w2[2]*w2[2];
    if ( dd1 > dd2 ) {
      dd1 = 1.0 / sqrt(dd1);
      v[0][0] = w1[0]*dd1;
      v[0][1] = w1[1]*dd1;
      v[0][2] = w1[2]*dd1;
    }
    else {
      dd2 = 1.0 / sqrt(dd2);
      v[0][0] = w2[0]*dd2;
      v[0][1] = w2[1]*dd2;
      v[0][2] = w2[2]*dd2;
    }

    /* 3rd vector orthogonal */
    v[1][0] = v[2][1]*v[0][2] - v[2][2]*v[0][1];
    v[1][1] = v[2][2]*v[0][0] - v[2][0]*v[0][2];
    v[1][2] = v[2][0]*v[0][1] - v[2][1]*v[0][0];
    dd1 = v[1][0]*v[1][0] + v[1][1]*v[1][1] + v[1][2]*v[1][2];
    dd1 = 1.0 / sqrt(dd1);
    v[1][0] *= dd1;
    v[1][1] *= dd1;
    v[1][2] *= dd1;
  }

  l[0] *= maxm;
  l[1] *= maxm;
  l[2] *= maxm;

  return(n);
}


/* 2d eigen value+vector extraction */
int eigen_2d(double *m,double *l,double vp[2][2]) {
  double   d,d1,d2;

  /* init */
  l[0] = l[1] = 0.0;
  vp[0][0] = 1.0;  vp[0][1] = 0.0;
  vp[1][0] = 0.0;  vp[1][1] = 1.0;

  /* eigenvalues of jacobian */
  if ( fabs(m[0]) > EV_EPSD ) {
    d = 0.5 * sqrt((m[0]-m[2])*(m[0]-m[2]) + 4*m[1]*m[1]);

    l[0] = 0.5 * (m[0]+m[2]) + d;
    l[1] = 0.5 * (m[0]+m[2]) - d;
    if ( fabs(l[0]) < EV_EPSD || fabs(l[1]) < EV_EPSD )  return(1);

    vp[0][0] = m[1];  vp[0][1] = l[0]-m[0];
    vp[1][0] = m[1];  vp[1][1] = l[1]-m[0];
  
    /* normalize */
    d1 = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);  
    d2 = sqrt(vp[1][0]*vp[1][0] + vp[1][1]*vp[1][1]);
    if ( d1 > EV_EPSD ) {
      vp[0][0] /= d1;
      vp[0][1] /= d1;
    }
    if ( d2 > EV_EPSD ) {
      vp[1][0] /= d2;
      vp[1][1] /= d2;
    }
  }

  return(1);
}
