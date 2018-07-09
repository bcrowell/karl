"""
A wrapper around the RK routines to automatically switch coordinate systems when
appropriate, and to accurately calculate the termination of incomplete geodesics.
"""

#include "language.h"
#include "util.h"
#include "math.h"
#include "io_util.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

#if "LANG" eq "python"
import ctypes,numpy
import c_libs
from c_libs import c_double_p
#endif

from io_util import fl

import schwarzschild,kruskal,angular,transform,runge_kutta

def trajectory_schwarzschild(spacetime,chart,x0,v0,opt,sigma):
  """
  A specialized wrapper for trajectory_simple(), optimized for speed and numerical precision in the
  Schwarzschild spacetime, for causal, future-oriented world-lines.
  Inputs and outputs are the same as for runge_kutta.trajectory_simple(), plus the additional input sigma,
  and opt['tol'], which is a desired error tolerance.

  Values of spacetime and chart are defined in spacetimes.h.
  x0 and v0 = initial position and velocity, expressed in the chosen chart
  opt = same hash of options defined in runge_kutta
  sigma = see docs, distinguishes regions I and II from III and IV, if this ambiguity isn't resolvable
             in the chosen coordinate chart
  returns [err,x,v,a,final_lambda,info,sigma]
  """
  ndim = 5
  ndim2 = 2*ndim
  # The following are needed only by rddot().
#if "LANG" eq "python"
  acc = numpy.zeros((ndim,))
  acc = acc.astype(numpy.float64)
  pt = numpy.zeros((ndim2,))
  pt = pt.astype(numpy.float64)
  acc_p = acc.ctypes.data_as(c_double_p)
  pt_p = pt.ctypes.data_as(c_double_p)
#else
  acc = EMPTY1DIM(ndim)
  pt = EMPTY1DIM(ndim2)
  acc_p = NONE
  pt_p = NONE
#endif
  if HAS_KEY(opt,'triggers'):
    user_triggers = CLONE_ARRAY_OF_FLOATS2DIM(opt['triggers'])
  else:
    user_triggers = []
  lambda0 = 0.0
  if HAS_KEY(opt,'lambda0'):
    lambda0 = opt['lambda0']
  real_lambda_max = opt['lambda_max'] # as opposed to the shorter chunks we do with fixed step size
  tol = opt['tol']
  x = CLONE_ARRAY_OF_FLOATS(x0)
  v = CLONE_ARRAY_OF_FLOATS(v0)
  user_chart = chart
  while True: #js while (true)
    triggers = CLONE_ARRAY_OF_FLOATS2DIM(user_triggers)
    # Find out our initial Schwarzschild r:
    xs = transform.transform_point(x,spacetime,chart,CH_SCH)
    r = xs[1]
    # Pick coordinates that are well adapted to the region we're currently in.
    if r>1.1:
      optimal_chart = CH_SCH
      APPEND_TO_ARRAY(triggers,([-1.0,1,1.05,0.3])) # trigger on r<1.05, nearing horizon
    else:
      if r>0.9:
        optimal_chart = CH_AKS
        APPEND_TO_ARRAY(triggers,([ 1.0,1, 0.0,0.3])) # trigger on b>0, entering region II
        APPEND_TO_ARRAY(triggers,([-1.0,1, 0.0,0.3])) # trigger on b<0, leaving horizon for region I
      else:
        if r>0.5:
          optimal_chart = CH_SCH
          APPEND_TO_ARRAY(triggers,([-1.0,1,0.35,0.3])) # trigger on r<0.05, nearing singularity
          APPEND_TO_ARRAY(triggers,([ 1.0,1,0.95,0.3]))
          # ... only relevant for a spacelike world-line
        else:
          optimal_chart = CH_KEP
          APPEND_TO_ARRAY(triggers,([1.0,1, 0.2,0.3]))
          # ... only relevant for a spacelike world-line
          # It's not useful to try to make a trigger that prevents hitting the singularity. Triggers are
          # too crude for that purpose, don't work reliably because the coordinate velocities diverge.
    if chart!=optimal_chart:
      x2 = transform.transform_point(x,spacetime,chart,optimal_chart)
      v2 = transform.transform_vector(v,x,spacetime,chart,optimal_chart)
      x = x2
      v = v2
    chart = optimal_chart
    opt['triggers'] = triggers
    # Use heuristics to pick a step size:
    # fixme -- the following only really applies at small r
    r_stuff = runge_kutta.r_stuff(spacetime,chart,x,v,acc,pt,acc_p,pt_p)
    r,rdot,rddot,p,lam_left = r_stuff
    # fixme -- sanity checks on lam_left
    # fixme: don't hardcode parameters here
    alpha = 0.5
    k = 0.3
    r_min = 1.0e-16
    # time to quit?
    if r<r_min:
      #PRINT("quitting because r<r_min")
      final_lambda = final_lambda+lam_left
      return final_helper(RK_INCOMPLETE,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,user_chart)
    # set step size
    dlambda = k*tol**0.25*lam_left**alpha
    n=100 # Try to do enough steps with fixed step size to avoid excessive overhead.
    safety = 0.7 # margin of safety so that we never hit singularity
    if n*dlambda>safety*lam_left:
      # ...fixme -- can this be improved?
      n=safety*lam_left/dlambda
      if n<1:
        n=1
        dlambda = safety*lam_left
    opt['dlambda'] = dlambda
    opt['lambda_max'] = lambda0+n*dlambda
#if 1
    PRINT("r=",io_util.fl(r),", r'=",io_util.fl(rdot),", r''=",io_util.fl(rddot),\
            ", lam=",io_util.fl(lambda0),\
            ", dlam=",io_util.fl(dlambda),", lam_max=",io_util.fl(opt['lambda_max']))
#endif

    result = runge_kutta.trajectory_simple(spacetime,chart,x,v,opt)
    err,final_x,final_v,final_a,final_lambda,info = result
#if 0
    PRINT("final_lambda=",final_lambda,", final_x=",io_util.vector_to_str_n_decimals(final_x,16))
#endif
    if chart==CH_AKS:
      sigma = kruskal.sigma(x[0],x[1]) # e.g., could have moved from III to II
    if final_lambda>=real_lambda_max:
      #PRINT("quitting because final_lambda>=real_lambda_max")
      return final_helper(err,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,user_chart)

    x = final_x
    v = final_v
    # Find final Schwarzschild r:
    xs = transform.transform_point(x,spacetime,chart,CH_SCH)
    r = xs[1]
    # Check for incomplete geodesic:
    if r<EPS or opt['dlambda']<EPS:
      r_stuff = runge_kutta.r_stuff(spacetime,chart,x,v,acc,pt,acc_p,pt_p)
      r,rdot,rddot,p,lam_left = r_stuff
      #PRINT("quitting because r<10*EPS or opt['dlambda']<10*EPS")
      return final_helper(RK_INCOMPLETE,final_x,final_v,final_a,final_lambda+lam_left,\
                          {'message':'incomplete geodesic'},\
                          sigma,spacetime,chart,user_chart)
    opt['triggers'] = user_triggers
    lambda0 = final_lambda
    opt['lambda0'] = final_lambda

def final_helper(err,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,desired_chart):
  x = transform.transform_point(final_x,spacetime,chart,desired_chart)
  v = transform.transform_vector(final_v,final_x,spacetime,chart,desired_chart)
  a = transform.transform_vector(final_a,final_x,spacetime,chart,desired_chart)
  return [err,x,v,a,final_lambda,info,sigma]

  
