import schwarzschild,util

import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

from sph_point import SphPoint
from sph_vector import SphVector

# x = starting point, a SphPoint object
# v = starting tangent vector (need not be normalized, can be null or spacelike)
# lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
#                  we hit a singularity)
# dlambda = step size
# returns
#   [if_error,error_message,final_x,final_v]
def sph_geodesic_rk(x,v,lambda_max,dlambda):
  spacetime = x.spacetime # SphPoint.SCHWARZSCHILD or SphPoint.CHARGED
  ok = False
  if spacetime==SphPoint.SCHWARZSCHILD: ok = True
  if spacetime==SphPoint.CHARGED: return [True,"CHARGED not implemented"]
  if not ok: return [True,"spacetime not implemented"]
  n = math.ceil(lambda_max/dlambda)
  for iter in range(0,n):
    if spacetime==SphPoint.SCHWARZSCHILD:
      theta = x.theta
      ch = schwarzschild.sch_christoffel_sch(x.t,x.r,sin(theta),cos(theta))
    a = [0.0,0.0,0.0,0.0] # second derivative of x^i with respect to lambda
    for i in range(0, 4):
      for j in range(0, 4):
        for k in range(0, 4):
          a[i] = a[i] + ch[j][k][i]*v.comp[j]*v.comp[k]
    for i in range(0, 4):
      v.comp[i] += dlambda*a[i]
    x.add(dlambda,v)
  return [False,"",x,v]
