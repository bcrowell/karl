#include "spacetimes_c.h"

#define NDIM 5
void sum_christoffel_sch_sch(double *,double *);

#define APPLY(a,p,i,j,k,ch) (a[i] -= ch*p[NDIM+j]*p[NDIM+k])

void apply_christoffel(int spacetime,int chart,double *p,double *a,double dlambda) {
  int m;
  for (m=0; m<NDIM; m++) {a[m]=0.0;}
  if ((spacetime | chart)==(SP_SCH | CH_SCH)) {
    sum_christoffel_sch_sch(p,a);
  }
  for (m=0; m<NDIM; m++) {a[m]=a[m]*dlambda;}
}

// For the Schwarzschild spacetime, described in 5-dimensional Schwarzschild coordinates,
// add up the Christoffel symbols as they appear in the geodesic equation.
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
  APPLY(a,p,  0,1,0,z);
  APPLY(a,p,  1,0,0,z);
  // These terms are analogous to Gamma^r_theta_theta=-c and Gamma^theta_r_theta=1/r in sch4:
  for (m=2; m<=4; m++) {
    APPLY(a,p,  m,m,1,-c);
    APPLY(a,p,  m,1,m,inv_r);
    APPLY(a,p,  1,m,m,inv_r);
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

