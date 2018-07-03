#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spacetimes_c.h"

#define NDIM 5
#define APPLY(a,p,i,j,k,ch) (a[k] -= (ch)*(p)[NDIM+(i)]*(p)[NDIM+(j)])
// ... third index is contravariant, as in ctensor

double veberic_lambert_w(double);
double lambert_w_of_exp(double);
void aux(double *,double *,double *,double,double);
void sum_christoffel_sch_sch(double *,double *);
void sum_christoffel_sch_aks(double *,double *);
void sum_christoffel_sch_kep(double *,double *);

// Given a position and velocity (both packed into the 10 elements of p),
// find the dv based on the Christoffel symbols and dlambda.
void apply_christoffel(int spacetime,int chart,double *p,double *a,double dlambda) {
  int m,ok;
  ok = 0;
  for (m=0; m<NDIM; m++) {a[m]=0.0;}
  if ((spacetime | chart)==(SP_SCH | CH_SCH)) {
    sum_christoffel_sch_sch(p,a);
    ok = 1;
  }
  if ((spacetime | chart)==(SP_SCH | CH_AKS)) {
    sum_christoffel_sch_aks(p,a);
    ok = 1;
  }
  if ((spacetime | chart)==(SP_SCH | CH_KEP)) {
    sum_christoffel_sch_kep(p,a);
    ok = 1;
  }
  if (!ok) {
    fprintf(stderr,"Unrecognized combination of spacetime=%d, chart=%d in apply_christoffel.c.\n",
                spacetime,chart);
    exit(-1);
  }
  for (m=0; m<NDIM; m++) {a[m]=a[m]*dlambda;}
}

// For the Schwarzschild spacetime, described in 5-dimensional Schwarzschild coordinates,
// add up the Christoffel symbols as they appear in the geodesic equation.
// The input array a[] needs to be zeroed before calling.
void sum_christoffel_sch_sch(double *p,double *a) {
  int m,n;
  double z,r,inv_r,r2,c;
  r = p[1];
  r2 = r*r;
  inv_r = 1.0/r;
  c = r-1.0; // = r-2GM in Carroll's notation
  APPLY(a,p,  0,0,1,0.5*c/(r2*r));
  z = 0.5/(r*c);
  APPLY(a,p,  1,1,1,-z);
  APPLY(a,p,  0,1,0,2.0*z); // factor of 2 handles 100 as well as 010
  // These terms are analogous to Gamma^r_theta_theta=-c and Gamma^theta_r_theta=1/r in sch4:
  for (m=2; m<=4; m++) {
    APPLY(a,p,  m,m,1,-c);
    APPLY(a,p,  m,1,m,2.0*inv_r); // factor of 2 handles 1mm as well as m1m
  }
  // Fictitious centripetal terms:
  // For efficiency, we assume xi^2=i^2+j^2+k^2=1.
  for (m=2; m<=4; m++) { // upper index
    z = p[m];
    for (n=2; n<=4; n++) { // lower indices
      APPLY(a,p,  n,n,m,z);
    }
  }
}

// For the Schwarzschild spacetime, described in 5-dimensional "Keplerian" coordinates,
// add up the Christoffel symbols as they appear in the geodesic equation.
// The input array a[] needs to be zeroed before calling.
void sum_christoffel_sch_kep(double *p,double *a) {
  int m,n;
  double z,u,c,c2,u2,u13,u23,u43,u53;
  u = p[1];
  u2 = u*u;
  u13 = pow(u,(1.0/3.0)); // u^(1/3)
  u23 = u13*u13; // u^(2/3)
  u43 = u23*u23; // u^(4/3)
  u53 = u43*u13; // u^(5/3)
  APPLY(a,p,  0,0,1,0.75*(1/u-1/u53)); // ^u _t t
  z = -(1.0/3.0)*u13/(u43-u2);
  APPLY(a,p,  0,1,0,2.0*z); // factor of 2 handles 100 as well as 010
  APPLY(a,p,  1,1,1,(1.0/3.0)*(1/u13)*((1.0-2.0*u23+u43)/(1-3.0*u23+3.0*u43-u2))); // ^u _u u
  // These terms are analogous to Gamma^r_theta_theta=-c and Gamma^theta_r_theta=1/r in sch4:
  c = (1.5)*(u13-u);
  c2 = (2.0/3.0)/u;
  for (m=2; m<=4; m++) {
    APPLY(a,p,  m,m,1,c);
    APPLY(a,p,  m,1,m,2.0*c2); // factor of 2 handles 1mm as well as m1m
  }
  // Fictitious centripetal terms:
  // For efficiency, we assume xi^2=i^2+j^2+k^2=1.
  for (m=2; m<=4; m++) { // upper index
    z = p[m];
    for (n=2; n<=4; n++) { // lower indices
      APPLY(a,p,  n,n,m,z);
    }
  }
}

// For the Schwarzschild spacetime, described in 5-dimensional arcsinh Kruskal coordinates,
// add up the Christoffel symbols as they appear in the geodesic equation.
// This mostly duplicates python code in kruskal.pp. It doesn't check as carefully for
// underflows in functions like exp() and tanh().
// The input array a[] needs to be zeroed before calling.
void sum_christoffel_sch_aks(double *p,double *a) {
  double aa,bb,t,r,mu,tanha,tanhb,q,q2,z;
  int m,n;
  aa=p[0];
  bb=p[1];
  aux(&t,&r,&mu,aa,bb);
  tanha = tanh(aa);
  tanhb = tanh(bb);
  q = 0.5*mu*(1.0+1.0/r);
  q2 = -0.5*mu/r;
  APPLY(a,p,  0,0,0,tanha+q*tanhb);
  APPLY(a,p,  1,1,1,tanhb+q*tanha);
  z = 2.0*q2*tanhb; // factor of 2 handles symmetry in the following
  APPLY(a,p,  0,2,2,z);
  APPLY(a,p,  0,3,3,z);
  APPLY(a,p,  0,4,4,z);
  z = 2.0*q2*tanha; // factor of 2 handles symmetry in the following
  APPLY(a,p,  1,2,2,z);
  APPLY(a,p,  1,3,3,z);
  APPLY(a,p,  1,4,4,z);
  z = -0.5*r*tanha;
  APPLY(a,p,  2,2,0,z);
  APPLY(a,p,  3,3,0,z);
  APPLY(a,p,  4,4,0,z);
  z = -0.5*r*tanhb;
  APPLY(a,p,  2,2,1,z);
  APPLY(a,p,  3,3,1,z);
  APPLY(a,p,  4,4,1,z);
  // Fictitious centripetal terms:
  // For efficiency, we assume xi^2=i^2+j^2+k^2=1.
  for (m=2; m<=4; m++) { // upper index
    z = p[m];
    for (n=2; n<=4; n++) { // lower indices
      APPLY(a,p,  n,n,m,z);
    }
  }
}

// A C implementation of kruskal.aux. See comments there.
// Inputs are Kruskal arcsinh coords (a,b), outputs are t,r,mu.
// Unlike the python implementation, this one doesn't check for
// coordinate values past the singularity, and doesn't try hard
// to avoid underflows in exponentials.
void aux(double *t,double *r,double *mu,double a,double b) {
  double e2a,e2b,f,u;
  if (a<0.0) {
    // Flip to positive a in order to simplify some later computations.
    aux(t,r,mu,-a,-b);
    return;
  }
  // From now on, we know a>=0.
  e2a = exp(-2.0*a);
  e2b = exp(-2.0*fabs(b));
  // From now on, we know we're in region I or II, a>0.
  f = (1.0-e2a)*(1.0-e2b);
  u = a+fabs(b)+log(f/4.0)-1.0;
  if (b<0.0) {
    // region I
    *r = 1.0+lambert_w_of_exp(u);
  }
  else {
    // region II
    *r = 1.0+veberic_lambert_w(-exp(u));
  }
  // Compute mu:
  *mu = (1.0+e2a)*(1.0+e2b)*(1/(2*M_E*(*r)))*exp(a+fabs(b)-((*r)-1.0));
  // Compute t:
  if (a!=0.0 && b!=0.0) {
    *t = a-fabs(b)+log((1.0-e2a)/(1.0-e2b));
  }
}

