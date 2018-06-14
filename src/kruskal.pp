"""
Low-level routines to compute things for points in the Schwarzschild
spacetime, in Kruskal-Szekeres coordinates, rescaled by an arcsinh,
with the angular space embedded on the unit 2-sphere in a fictitious
additional dimension.  The units are such that the Schwarzschild
radius is 1, i.e., the mass is 1/2. 

Transformations between arcsinh-Kruskal coordinates and Schwarzschild
coordinates, and the associated jacobians, are in transform.pp, except
that transformation from (a,b) to (t,r) are here, in kruskal.aux.

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.) 
"""

#include "util.h"
#include "math.h"
#include "precision.h"

import math_util
from math_util import safe_exp

def metric_ks4(p):
  """
  For the Schwarzschild spacetime, compute the metric.
  Metric is in lower-index form, in +--- signature, with coordinates (a,b,theta,phi).
  """
  g = [[0 for i in range(4)] for j in range(4)]
  t,r,mu = aux(p[0],p[1])
  g[0][1] = mu
  g[1][0] = mu
  r2 = r*r
  g[0][0] = a
  g[1][1] = -1.0/a
  g[2][2] = -r2
  g[3][3] = -r2*sin_theta*sin_theta
  return g


def aux(a,b):
  """
  Compute Schwarzschild t and r, and metric element mu, for a point given in rescaled Kruskal coordinates (a,b).
  Return [t,r,mu].

  Although r is an analytic function of the coordinates throughout the entire
  maximal extension of the spacetime, numerical issues force us to do some case-splitting.
  The quantity mu is the coefficient in the line element, ds^2 = 2 mu da db -...,
  i.e., it's the matrix element g_ab of the metric. Returns t=None for points on a horizon.
  Returns None for all three variables if a and b lie beyond the singularity.
  """
  if a<=0: # Flip to positive a in order to simplify some later computations.
    return aux(-a,-b)
  # From now on, we know a>=0.
  e2a = safe_exp(-2*a)
  e2b = safe_exp(-2*abs(b))
  if b>=0.0: # Can I simplify the logic by eliminating this?
    # interior (region II) or horizon of II
    r,mu = small_r_mu(a,b)
  else:
    # exterior, region I
    # From now on, we know a>0 and b<0.
    # Here, if r is even moderately large the coordinates v and w can be too big to store in
    # floating point, so we use a different algorithm.
    f = (1.0-e2a)*(1.0-e2b)
    u = a-b+log(f/4.0)-1.0
    r = 1.0+lambert_w_of_exp(u)
    # Compute mu:
    mu = (1.0+e2a)*(1.0+e2b)*(1/(2*math.e*r))*exp(a-b-(r-1))
  # Compute t:
  if a!=0 and b!=0 and r is not None:
    t = a-abs(b)+log((1-e2a)/(1-e2b))
  else:
    t = None
  return [t,r,mu]

def small_r_mu(a,b):
  """
  Returns [r,mu]. Works only at small r and t, produces overflows elsewhere. Returns [None,None] in
  regions beyond the singularity.
  """
  v=sinh(a); w=sinh(b)
  rho = -v*w
  if rho<-1.0: return [None,None]
  r = 1+lambert_w(rho/math.e)
  mu = 2.0*((r-1)/(r*rho))*cosh(a)*cosh(b)
  return [r,mu]

def christoffel(p):
  """
  For the Schwarzschild spacetime, compute the Christoffel symbols, in coordinates (a,b,i,j,k).

  The order of indices is as in ctensor:
     symmetric on 1st 2 indices
     contravariant on final index
  See maxima/schwarzschild5.mac. In addition, we put in a fictitious centripetal term that
  keeps us constrained to the unit sphere i^2+j^2+k^2=1.
  """
  ch = [[[0 for i in range(5)] for j in range(5)] for k in range(5)]
  r,mu = r(p[0],p[1])
  # In the following, need to insert the appropriate Christoffel symbols, including the
  # terms analogous to Gamma^r_theta_theta.
  raise RuntimeError('FIXME') # not yet implemented
  # Fictitious centripetal terms:
  i=p[2] ; j=p[3]; k=p[4]
  xi2 = i*i+j*j+k*k; # should normally be very close to 1
  for m in range(2,5): # upper index
    z = p[m]
    for n in range(2,5): # lower indices
      ch[n][n][m] = z/xi2
  return ch

def christoffel_sch4(t,r,sin_theta,cos_theta):
  """
  For the Schwarzschild spacetime, compute the Christoffel symbols, in Schwarzschild coordinates (t,r,theta,phi).

  The order of indices is as in ctensor:
     symmetric on 1st 2 indices
     contravariant on final index
  This will give a division by zero exception if used at the coordinate singularities at theta=0 and pi
  and r=1.
  The equations are from http://ned.ipac.caltech.edu/level5/March01/Carroll3/Carroll7.html ,
  equation 7.33.
  """
  ch = [[[0 for i in range(4)] for j in range(4)] for k in range(4)]
  r2 = r*r
  c = r-1.0 # = r-2GM in Carroll's notation
  ch[0][0][1] = 0.5*c/(r2*r)
  z = 0.5/(r*c)
  ch[1][1][1] = -z
  ch[0][1][0] = z
  ch[1][0][0] = z
  z = 1.0/r
  ch[1][2][2] = z
  ch[2][1][2] = z
  ch[2][2][1] = -c
  ch[1][3][3] = z
  ch[3][1][3] = z
  s2 = sin_theta*sin_theta
  ch[3][3][1] = -c*s2
  ch[3][3][2] = -sin_theta*cos_theta
  cot = cos_theta/sin_theta
  ch[2][3][3] = cot
  ch[3][2][3] = cot
  return ch

