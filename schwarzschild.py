#!/usr/bin/python3

# Low-level routines to compute things for the Schwarzschild spacetime.

# Documentation for the math is in the file doc.tex, which can be 
# compiled to pdf format by doing a "make doc." (Comments in the code do
# not document the math or the definitions of the variables.)


import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

from util import lambert_w

# Define coordinate charts for different "eras," i.e., ranges of Schwarzschild t.
# Representation is (t,V,W), where t represents a time translation to be applied
# to the point (V,W), and we try to keep (V,W) such that it is fairly close to t=0.
# This prevents (V,W) from getting too big to represent in floating point. (They grow 
# exponentially with |t|.)
# This routine checks whether we've drifted out of the canonical chart, and if we have,
# fixes it.
def ks_era_in_range(ks):
  t = ks[0]
  v = ks[1]
  w = ks[2]
  tx = ks_tx(v,w)
  if abs(tx)<5.0:
    return [False,ks] # no change needed, already in canonical chart
  else:
    return [True,force_ks_era(ks,tx,ks_to_region(v,w))]

# Force (t,V,W) into a canonical form where (V,W) corresponds to a Schwarzschild t=0,
# i.e., the t information is all in the t parameter.
# See ks_tx() for definition of tx. See ks_era_jacobian() for the corresponding
# Jacobian matrix.
def force_ks_era(ks,tx,region):
  t = ks[0]
  v = ks[1]
  w = ks[2]
  t = t + 2.0*arctanh(tx)
  root_rho = sqrt(abs(v*w)) # sqrt of absolute value of rho
  v=root_rho*sign(v)
  w=root_rho*sign(w)
  return [t,v,w]

# This can cause numerical overflows, should typically only be attempted for output.
def ks_to_zero_era(t,v,w):
  tr = ks_to_sch(v,w)
  return sch_to_ks(tr[0]+t,tr[1],ks_to_sigma(v,w))

# calculate T/X if exterior or X/T if interior:
def ks_tx(v,w):
  region = ks_to_region(v,w)
  is_exterior = ks_is_exterior(region)
  if is_exterior:
    return (v+w)/(v-w)
  else:
    return (v-w)/(v+w)

def sch_to_sigma(r):
  if r>=1.0:
    return 1
  else:
    return -1

def ks_to_region(v,w):
  if v>=0 and w<=0: return 1 # exterior
  if v>=0 and w>=0: return 2 # interior
  if v<=0 and w>=0: return 3 # exterior, other copy of Minkowski space
  return 4 # white hole

def ks_is_exterior(region):
  return (region==1 or region==3)

def ks_to_sigma(v,w):
  return 1.0 if v>0 else -1.0  

# Take null Kruskal-Szekeres coordinates (V,W,theta,phi) and convert them 
# to Schwarzschild coordinates (t,r,theta,phi).
# The V,W coordinates are equivalent to Hawking and Ellis's (v'/sqrt2,w'/sqrt2).
# This will inevitably misbehave at the horizon, since KS coordinate are eliminating
# a coordinate singularity.
# This is my own way of expressing this transformation, may contain mistakes, needs to be tested.
def ks_to_sch(v,w):
  rho = -v*w
  l = lambert_w(rho/math.e)
  r = 1+l
  ks_x = (v-w)/2.0
  ks_t = (v+w)/2.0
  if rho>0.0:
    t = 2.0*arctanh(ks_t/ks_x) # exterior
  else:
    t = 2.0*arctanh(ks_x/ks_t) # interior
  return [t,r]

# This is guaranteed to misbehave at the horizon.
# sigma=+1 for regions I and II, -1 for III and IV
def sch_to_ks(t,r,sigma):
  t2 = 0.5*t
  if r>1.0:
    region=1 if sigma>1.0 else 3
    sc = sinh(t2)
    cs = cosh(t2)
  else:
    region=2 if sigma>1.0 else 4
    sc = cosh(t2)
    cs = sinh(t2)
  # My formulation based on MTW p. 827:
  a = sigma*sqrt(abs(r-1.0))*exp(r/2.0)
  ks_t = a*sc
  ks_x = a*cs
  v = ks_t+ks_x
  w = ks_t-ks_x
  return [v,w]

# For the Schwarzschild spacetime, compute the Christoffel symbols, in Schwarzschild coordinates
# (t,r,theta,phi).
# The order of indices is as in ctensor:
#    symmetric on 1st 2 indices
#    contravariant on final index
# This will give a division by zero exception if used at the coordinate singularities at theta=0 and pi
# and r=1.
# The equations are from http://ned.ipac.caltech.edu/level5/March01/Carroll3/Carroll7.html ,
# equation 7.33.
def sch_christoffel_sch(t,r,sin_theta,cos_theta):
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

# For the Schwarzschild spacetime, compute the Christoffel symbols, in Kruskal-Szekeres null coordinates
# (V,W,theta,phi).
# The order of indices is as in ctensor:
#    symmetric on 1st 2 indices
#    contravariant on final index
# See docs for readable mathematical versions of these.
# This will give a division by zero exception if used at the coordinate singularities at theta=0 and pi.
def sch_christoffel_ks(v,w,sin_theta,cos_theta,r,b):
  ch = [[[0 for i in range(4)] for j in range(4)] for k in range(4)]
  qb = 0.25*b   # = 1/(re^r)
  qbr = qb/r    # = 1/(r^2e^r)
  c = qb*(1.0+1.0/r)
  ch[0][0][0] = c*w        # _VV ^V
  ch[1][1][1] = c*v        # _WW ^W
  z = -w*qbr
  ch[0][2][2] = z          # _Vtheta ^theta
  ch[2][0][2] = z
  ch[0][3][3] = z          # _Vphi   ^phi
  ch[3][0][3] = z
  z = -v*qbr
  ch[1][2][2] = z          # _Wtheta ^theta
  ch[2][1][2] = z
  ch[1][3][3] = z          # _Wphi   ^phi
  ch[3][1][3] = z
  eighth_r = 0.125*r
  ch[2][2][0] = -v*eighth_r  # _thetatheta ^V
  ch[2][2][1] = -w*eighth_r  # _thetatheta ^W
  sin2_theta = sin_theta*sin_theta
  ch[3][3][0] = -v*eighth_r*sin2_theta  # _phiphi ^V
  ch[3][3][1] = -w*eighth_r*sin2_theta  # _phiphi ^W
  ch[3][3][2] = -sin_theta*cos_theta
  z = cos_theta/sin_theta # will give an exception at theta=0, pi
  ch[2][3][3] = z          # _thetaphi ^phi
  ch[3][2][3] = z
  return ch

# For testing purposes, compute Christoffel symbols for the Schwarzschild spacetime in Kruskal
# coordinates, using the (nearly) raw output from Maxima. The code below was produced by executing
# christoffel.mac and then running clean_up_christoffel.rb.
def sch_christoffel_ks_raw_maxima(v,w,theta,r):
  ch = [[[0 for i in range(4)] for j in range(4)] for k in range(4)]
  ch[0][0][0] = (-w/((r-1)*exp((r-1)+1)+2*(r-1)**2*exp((r-1)+1)+(r-1)**3*exp((r-1)+1)))-(r-1)**2/(v+2*(r-1)*v+(r-1)**2*v)-(2*(r-1))/(v+2*(r-1)*v+(r-1)**2*v)-1/(v+2*(r-1)*v+(r-1)**2*v) 
  #   ... ^v _v v
  ch[0][2][2] = -w/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^theta _v theta
  ch[2][0][2] = -w/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^theta _theta v
  ch[0][3][3] = -w/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^phi _v phi
  ch[3][0][3] = -w/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^phi _phi v
  ch[1][1][1] = (-(r-1)**2/(w+2*(r-1)*w+(r-1)**2*w))-(2*(r-1))/(w+2*(r-1)*w+(r-1)**2*w)-1/(w+2*(r-1)*w+(r-1)**2*w)-v/((r-1)*exp((r-1)+1)+2*(r-1)**2*exp((r-1)+1)+(r-1)**3*exp((r-1)+1)) 
  #   ... ^w _w w
  ch[1][2][2] = -v/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^theta _w theta
  ch[2][1][2] = -v/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^theta _theta w
  ch[1][3][3] = -v/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^phi _w phi
  ch[3][1][3] = -v/(exp((r-1)+1)+2*(r-1)*exp((r-1)+1)+(r-1)**2*exp((r-1)+1)) 
  #   ... ^phi _phi w
  ch[2][2][0] = (exp((-(r-1))-1)*v**2*w)/(8*(r-1))+(exp((-(r-1))-1)*v**2*w)/8 
  #   ... ^v _theta theta
  ch[2][2][1] = (exp((-(r-1))-1)*v*w**2)/(8*(r-1))+(exp((-(r-1))-1)*v*w**2)/8 
  #   ... ^w _theta theta
  ch[2][3][3] = cos(theta)/sin(theta) 
  #   ... ^phi _theta phi
  ch[3][2][3] = cos(theta)/sin(theta) 
  #   ... ^phi _phi theta
  ch[3][3][0] = (exp((-(r-1))-1)*sin(theta)**2*v**2*w)/(8*(r-1))+(exp((-(r-1))-1)*sin(theta)**2*v**2*w)/8 
  #   ... ^v _phi phi
  ch[3][3][1] = (exp((-(r-1))-1)*sin(theta)**2*v*w**2)/(8*(r-1))+(exp((-(r-1))-1)*sin(theta)**2*v*w**2)/8 
  #   ... ^w _phi phi
  ch[3][3][2] = -cos(theta)*sin(theta)
  #   ... ^theta _phi phi
  return ch

# Compute the auxiliary quantities [rho,r,b] for a point given in Kruskal-Szekeres null coordinates.
# See docs for definitions of these quanties.
def sch_aux_ks(v,w):
  rho = -v*w
  r = 1+lambert_w(rho/math.e)
  b = 4.0*(r-1)/(r*rho)
  return [rho,r,b]

# For the Schwarzschild spacetime, compute the metric, in Kruskal-Szekeres null coordinates.
# Metric is in lower-index form, in +--- signature, with coordinates (V,W,theta,phi).
# See docs for the quantity B, which can be calculated using sch_aux_ks.
# The V,W coordinates are equivalent to Hawking and Ellis's (v'/sqrt2,w'/sqrt2).
def sch_metric_ks(v,w,sin_theta,r,b):
  g = [[0 for i in range(4)] for j in range(4)]
  half_b =  0.5*b  
  g[0][1] = half_b
  g[1][0] = half_b 
  r2 = r*r
  g[2][2] = -r2
  g[3][3] = -r2*sin_theta*sin_theta
  return g

# For the Schwarzschild spacetime, compute the metric, in Schwarzschild coordinates.
# Metric is in lower-index form, in +--- signature, with coordinates (t,r,theta,phi).
# The mass is assumed to be 1/2, so that r is in units of the Schwarzschild radius.
# Angles in radians.
# Metric taken from https://en.wikipedia.org/wiki/Schwarzschild_metric .
def sch_metric_sch(r,sin_theta):
  g = [[0 for i in range(4)] for j in range(4)]
  a = 1.0-1.0/r
  r2 = r*r
  g[0][0] = a
  g[1][1] = -1.0/a
  g[2][2] = -r2
  g[3][3] = -r2*sin_theta*sin_theta
  return g

# Find the Jacobian matrix for the Kruskal V,W coordinates as functions of Schwarzschild t,r.
# The following was based on code generated by maxima, see Makefile in maxima directory.
# The first index is 0 for V, 1 for W. The second index is 0 for t, 1 for r.
# These are to be used in transitioning a vector from Sch. chart to Kruskal.
# Some of these can be simplified quite a bit if you know the VW coordinates, e.g., in region 1,
# jacobian[0][0]=V/2 and jacobian[1][0]=W/2. Probably not worthwhile to fiddle with optimizing
# that, because typically we only use this when transitioning from one chart to another, which
# happens O(1) times in a calculation.
def sch_ks_jacobian(region,r,t):
  er = exp(r/2.0)
  c = cosh(t/2.0)
  s = sinh(t/2.0)
  jacobian = [[0 for i in range(2)] for j in range(2)]
  if region==1:
    sr1 = sqrt(r-1.0)
    jacobian[0][0] =  (sr1*er*s)/2+(sr1*er*c)/2 
    jacobian[0][1] =  (sr1*er*s)/2+(er*s)/(2*sr1)+(sr1*er*c)/2+(er*c)/(2*sr1) 
    jacobian[1][0] =  (sr1*er*c)/2-(sr1*er*s)/2 
    jacobian[1][1] =  (sr1*er*s)/2+(er*s)/(2*sr1)-(sr1*er*c)/2-(er*c)/(2*sr1) 
  if region==2:
    s1r = sqrt(1.0-r)
    jacobian[0][0] =  (s1r*er*s)/2+(s1r*er*c)/2 
    jacobian[0][1] =  (s1r*er*s)/2-(er*s)/(2*s1r)+(s1r*er*c)/2-(er*c)/(2*s1r) 
    jacobian[1][0] =  (s1r*er*s)/2-(s1r*er*c)/2 
    jacobian[1][1] =  (-(s1r*er*s)/2)+(er*s)/(2*s1r)+(s1r*er*c)/2-(er*c)/(2*s1r) 
  if region==3:
    sr1 = sqrt(r-1.0)
    jacobian[0][0] =  (sr1*er*s)/2+(sr1*er*c)/2 
    jacobian[0][1] =  (sr1*er*s)/2+(er*s)/(2*sr1)+(sr1*er*c)/2+(er*c)/(2*sr1) 
    jacobian[1][0] =  (sr1*er*c)/2-(sr1*er*s)/2 
    jacobian[1][1] =  (sr1*er*s)/2+(er*s)/(2*sr1)-(sr1*er*c)/2-(er*c)/(2*sr1) 
  if region==4:
    s1r = sqrt(1.0-r)
    jacobian[0][0] =  (s1r*er*s)/2+(s1r*er*c)/2 
    jacobian[0][1] =  (s1r*er*s)/2-(er*s)/(2*s1r)+(s1r*er*c)/2-(er*c)/(2*s1r) 
    jacobian[1][0] =  (s1r*er*s)/2-(s1r*er*c)/2 
    jacobian[1][1] =  (-(s1r*er*s)/2)+(er*s)/(2*s1r)+(s1r*er*c)/2-(er*c)/(2*s1r) 
  return jacobian

# Similar to sch_ks_jacobian.
# Find the Jacobian matrix for the Schwarzschild t,r coordinates as functions of Kruskal V,W.
# Slightly inefficient because we invert the 2x2 matrix,
# but not likely to affect performance because not called often.
# This can throw an error if called for a point on the horizon (r=1), where
# the Schwarzschild coordinates misbehave, and in general it's probably
# not going to be numerically accurate to use this near the horizon.
def ks_sch_jacobian(region,r,t):
  if r==1.0: raise RuntimeError('r=1 in ks_sch_jacobian')
  j = sch_ks_jacobian(region,r,t)
  a = j[0][0]
  d = j[1][1]
  b = j[1][0]
  c = j[0][1]
  det = a*d-b*c # should be nonzero because we checked above for r==1
  j[1][1] = a/det
  j[0][0] = d/det
  j[1][0] = -b/det
  j[0][1] = -c/det
  return j

# Calculate the jacobian matrix for the transformation of (t,V,W) to (0,V',W')
# implemented by force_ks_era(). The jacobian is 2x2, i.e., it doesn't describe
# what happens to t, and therefore it is always degenerate. The inputs v and
# w are the pre-transformation coordinates (V,W). The Jacobian has the form
#     W W
#     V V,
# multiplied by a scalar. The order of the indices is such that, e.g.,
#   dV'/dW = j_01.
def ks_era_jacobian(v,w):
  jacobian = [[0 for i in range(2)] for j in range(2)]
  rho = -v*w # is the same pre- and post-transformation
  # s=sign(rho)
  if rho>0.0:
    s = 1
  else:
    s = -1.0
  a = -s*0.5*(1/sqrt(abs(rho)))
  jacobian[0][0] = a*w
  jacobian[1][1] = a*v
  jacobian[1][0] = a*w
  jacobian[0][1] = a*v

