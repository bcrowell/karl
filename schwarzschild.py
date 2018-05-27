#!/usr/bin/python3

# Low-level routines to compute things for the Schwarzschild spacetime.

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
    return ks # no change needed, already in canonical chart
  else:
    return force_ks_era(ks,tx,ks_to_region(v,w))

# Force (t,V,W) into a canonical form where (V,W) corresponds to a Schwarzschild t=0,
# i.e., the t information is all in the t parameter.
def force_ks_era(ks,tx,region):
  t = ks[0]
  v = ks[1]
  w = ks[2]
  t = t + 2.0*arctanh(tx)
  is_exterior = ks_is_exterior(region)
  if is_exterior:
    root_rho = sqrt(-v*w)
  else:
    root_rho = sqrt(v*w)
  v=root_rho*sign(v)
  w=root_rho*sign(w)
  return [t,v,w]

def ks_tx(v,w):
  region = ks_to_region(v,w)
  is_exterior = ks_is_exterior(region)
  # calculate T/X if exterior or X/T if interior:
  if is_exterior:
    return (v+w)/(v-w)
  else:
    return (v-w)/(v+w)


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

# For the Schwarzschild spacetime, compute the metric, given Kruskal-Szekeres null coordinates.
# Metric is in lower-index form, in +--- signature, with coordinates (V,W,theta,phi).
# The V,W coordinates are equivalent to Hawking and Ellis's (v'/sqrt2,w'/sqrt2).
def sch_metric_ks(v,w,theta,phi):
  g = [[0 for i in range(4)] for j in range(4)]
  rho = -v*w
  l = lambert_w(rho/math.e)
  r = 1+l
  b = 4.0*l/((1+l)*rho)
  g[0][1] = 0.5*b  
  g[1][0] = 0.5*b  
  r2 = r**2
  g[2][2] = -r2
  g[3][3] = -r2*sin(theta)**2
  return g

# For the Schwarzschild spacetime, compute the metric, given Schwarzschild coordinates.
# Metric is in lower-index form, in +--- signature, with coordinates (t,r,theta,phi).
# The mass is assumed to be 1/2, so that r is in units of the Schwarzschild radius.
# Angles in radians.
# Metric taken from https://en.wikipedia.org/wiki/Schwarzschild_metric .
def sch_metric_sch(r,theta):
  g = [[0 for i in range(4)] for j in range(4)]
  a = 1.0-1.0/r
  r2 = r**2
  g[0][0] = a
  g[1][1] = -1.0/a
  g[2][2] = -r2
  g[3][3] = -r2*sin(theta)**2
  return g


