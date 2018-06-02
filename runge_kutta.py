import copy

import numpy
numpy.seterr(all='raise')
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

import schwarzschild,util

import io_util

# Calculate a geodesic using geodesic equation and 4th-order Runge-Kutta.
# Transitions between charts are made automatically as necessary.
# x = starting point, a SphPoint object or a member of some other class that implements a similar interface
# v = starting tangent vector (need not be normalized, can be null or spacelike); a SphVector object or
#     member of some other class that implements a similar interface
# lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
#                  we hit a singularity)
# dlambda = step size
# ndebug = 0, or, if nonzero, determines how often to print debugging output; e.g., if ndebug=100,
#           then we print debugging output after every ntrans*100 steps
# ndebug_inner -- similar to ndebug, but at a finer level of granularity; if it's 10, then
#            we print debugging info at every 10 steps
# ntrans -- after this many steps, we check whether to transition to a different chart; as long as
#           ntrans is at least 10, the time-consuming code for checking for a transition will only
#           run 1/10 of the time, and is therefore probably not a significant use of CPU time
# debug_transitions -- boolean, do we want a message printed whenever we make a transition?
# returns
#   [if_error,error_message,final_x,final_v,final_lambda]
def geodesic_rk(x,v,lambda_max,dlambda,ndebug,ndebug_inner,ntrans,debug_transitions):
  ok = False
  n = math.ceil(lambda_max/(dlambda*ntrans))
  if ndebug==0:
    steps_between_debugging=n*2 # debugging will never happen
  else:
    steps_between_debugging=n/ndebug
  debug_count = steps_between_debugging+1 # trigger it on the first iteration
  lam = 0.0
  x.make_safe()
  for iter in range(0,n):
    coords = x.get_raw_coords()
    do_debug = False
    if ndebug!=0 and (debug_count>=steps_between_debugging or iter==n-1):
      debug_count = 0
      do_debug = True
    if do_debug:
      z = x.absolute_schwarzschild()
      print("iter=",iter,", lambda=",("%5.3e" % lam),", chart=",x.chart,
            " coords=",io_util.vector_to_str(coords)," t=",("%3.1e" % z[0])," r=",("%3.1e" % z[1]))
    result = geodesic_rk_simple(x,v,lambda_max/n,dlambda,ndebug_inner)
    err = result[0]
    if err: return result
    #---- update the data
    x = result[2]
    v = result[3]
    lam = lam+result[4]
    result[4] = lam
    #---- Check for a transition to a new chart
    x.make_safe()
    debug_count += 1
  return result

# Calculate a geodesic using geodesic equation and 4th-order Runge-Kutta.
# The motion is assumed to stay within a single chart. If you want to do motion that may move from one
# chart to another, don't call this routine directly. Instead, call the fancier routine 
# geodesic_rk that acts as a wrapper for this one. Similar considerations apply if the
# motion may hit a singularity. Inputs are the same as the inputs to geodesic_rk of the same names.
# returns
#   [if_error,error_message,final_x,final_v,final_lambda]
def geodesic_rk_simple(x,v,lambda_max,dlambda,ndebug):
  ok = False
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
      z = x.absolute_schwarzschild()
      print("iter=",iter,", lambda=",("%5.3e" % lam),", chart=",x.chart,
            " coords=",io_util.vector_to_str(coords)," r=",("%5.3e" % z[1]))
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
      ch = x.get_christoffel(y)
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
