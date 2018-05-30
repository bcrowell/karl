import copy

import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

import schwarzschild,util

import io_util
from sph_point import SphPoint
from sph_vector import SphVector

# Calculate a geodesic using geodesic equation and 4th-order Runge-Kutta.
# x = starting point, a SphPoint object
# v = starting tangent vector (need not be normalized, can be null or spacelike)
# lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
#                  we hit a singularity)
# dlambda = step size
# ndebug = 0, or, if nonzero, gives approximate number of points at which to print debugging output
# returns
#   [if_error,error_message,final_x,final_v,final_lambda]
def sph_geodesic_rk(x,v,lambda_max,dlambda,ndebug):
  spacetime = x.spacetime # SphPoint.SCHWARZSCHILD or SphPoint.CHARGED
  ok = False
  if spacetime==SphPoint.SCHWARZSCHILD: ok = True
  if spacetime==SphPoint.CHARGED: return [True,"CHARGED not implemented"]
  if not ok: return [True,"spacetime not implemented"]
  n = math.ceil(lambda_max/dlambda)
  if ndebug==0:
    steps_between_debugging=n*2 # debugging will never happen
  else:
    steps_between_debugging=n/ndebug
  debug_count = steps_between_debugging+1 # trigger it on the first iteration
  lam = 0.0
  for iter in range(0,n):
    est = [[0 for i in range(8)] for step in range(4)]
            # four estimates of the changes in the independent variables for 4th-order Runge-Kutta 
            # reduce 2nd-order ODE to 8 coupled 1st-order ODEs
            # =k in the notation of most authors
    coords = x.get_raw_coords()
    do_debug = False
    if ndebug!=0 and (debug_count>=steps_between_debugging or iter==n-1):
      debug_count = 0
      do_debug = True
    if do_debug:
      print("iter=",iter,", lambda=",("%5.3e" % lam),", coords=",io_util.vector_to_str(coords))
    debug_count += 1
    y0 = [0 for i in range(8)]
    for i in range(0,4): y0[i]=copy.deepcopy(coords[i])
    for i in range(0,4): y0[i+4]=copy.deepcopy(v.comp[i])
    for step in range(0,4):
      if step==0: y=copy.deepcopy(y0)
      if step==1:
        for i in range(0,8):
          y[i] = y0[i]+0.5*est[0][i]
      if step==2:
        for i in range(0,8):
          y[i] = y0[i]+0.5*est[1][i]
      if step==3:
        for i in range(0,8):
          y[i] = y0[i]+est[2][i]
      if spacetime==SphPoint.SCHWARZSCHILD:
        theta = y[2]
        if x.chart==SphPoint.SCHWARZSCHILD:
          ch = schwarzschild.sch_christoffel_sch(y[0],y[1],sin(theta),cos(theta))
        else:
          return [True,"Kruskal chart not implemented"]
      # no else here because we tested at top to make sure it was Sch. spacetime
      for i in range(0,8): est[step][i]=0.0
      for i in range(0, 4):
        a = 0.0 # is essentially the acceleration
        for j in range(0, 4):
          for k in range(0, 4):
            a -= ch[j][k][i]*y[4+j]*y[4+k]
        est[step][i+4] = a*dlambda
      for i in range(0, 4):
        est[step][i] = y[4+i]*dlambda
    lam= lam+dlambda
    tot_est = [0 for i in range(8)]
    for i in range(0,8):
      tot_est[i] = (est[0][i]+2.0*est[1][i]+2.0*est[2][i]+est[3][i])/6.0
    for i in range(0, 4):
      v.comp[i] += tot_est[i+4]
    for i in range(0, 4):
      coords[i] += tot_est[i]    
    x.set_raw_coords(coords)
  return [False,"",x,v,lam]
