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

import schwarzschild,kruskal,angular,transform,runge_kutta,conserved

def trajectory_schwarzschild(spacetime,chart,x0,v0,opt):
  """
  A specialized wrapper for trajectory_simple(), optimized for speed and numerical precision in the
  Schwarzschild spacetime, for causal world-lines.
  Inputs and outputs are the same as for runge_kutta.trajectory_simple(), with some additional
  information packed in the options hash.

  Values of spacetime and chart are defined in spacetimes.h.
  x0 and v0 = initial position and velocity, expressed in the chosen chart
  opt = same hash of options defined in runge_kutta, plus some additional parameters:--
    required parameters:
      tol: desired error tolerance, notated delta in the docs
      sigma: see docs, distinguishes regions I and II from III and IV, if this ambiguity isn't resolvable
             in the chosen coordinate chart
      future_oriented: boolean; e.g., if we're in Schwarzschild coordinates, r<1, sigma==-1, and
             future_oriented is true, then we're rising away from the white-hole singularity of region III
    optional parameters:
      triggers: a set of triggers; if this option is supplied, allow_transitions must be false
      allow_transitions: boolean, should we allow transitions between coordinate systems?; default: true
  returns [err,x,v,a,final_lambda,info,sigma]
  """
  #---------------- Initialize some data. ------------------
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
  x = CLONE_ARRAY_OF_FLOATS(x0)
  v = CLONE_ARRAY_OF_FLOATS(v0)
  user_chart = chart
  #---------------- Read options. ------------------
  if HAS_KEY(opt,'allow_transitions'):
    allow_transitions = opt['allow_transitions']
  else:
    allow_transitions = TRUE
  if HAS_KEY(opt,'triggers'):
    if allow_transitions:
      THROW('allow_transitions should be false if there are user-supplied triggers')
    user_triggers = CLONE_ARRAY_OF_FLOATS2DIM(opt['triggers'])
  else:
    user_triggers = []
  lambda0 = 0.0
  if HAS_KEY(opt,'lambda0'):
    lambda0 = opt['lambda0']
  real_lambda_max = opt['lambda_max'] # as opposed to the shorter chunks we do with fixed step size
  tol = opt['tol']
  sigma = opt['sigma']
  future_oriented = opt['future_oriented']
  #---------------- Loop in which we check for coordinate transitions, then delegate the
  #                 real work to the basic RK routine. ------------------
  while True: #js while (true)
    triggers = CLONE_ARRAY_OF_FLOATS2DIM(user_triggers)
    r_stuff = runge_kutta.r_stuff(spacetime,chart,x,v,acc,pt,acc_p,pt_p)
    err,r,rdot,rddot,p,lam_left = r_stuff
    if err!=0:
      THROW("error in r_stuff; this probably means we hit the singularity, adaptive RK not working right")
    if allow_transitions:
      optimal_chart = chart_and_triggers(r,triggers,sigma,future_oriented)
      if chart!=optimal_chart:
        x2 = transform.transform_point(x,spacetime,chart,optimal_chart)
        v2 = transform.transform_vector(v,x,spacetime,chart,optimal_chart)
        x = x2
        v = v2
      chart = optimal_chart
    opt['triggers'] = triggers
    n,dlambda,terminate = choose_step_size(r,p,tol,lam_left,x,v,chart)
    if terminate:
      final_lambda = final_lambda+lam_left
      err = RK_INCOMPLETE
      BREAK
    opt['dlambda'] = dlambda
    opt['lambda_max'] = lambda0+n*dlambda
#if 0
    PRINT("r=",io_util.fl(r),", r'=",io_util.fl(rdot),", r''=",io_util.fl(rddot),\
            ", lam=",io_util.fl(lambda0),\
            ", dlam=",io_util.fl(dlambda),", lam_max=",io_util.fl(opt['lambda_max']))
#endif
    #---------------- Do n iterations of RK. ------------------
    result = runge_kutta.trajectory_simple(spacetime,chart,x,v,opt)
    err,final_x,final_v,final_a,final_lambda,info = result
#if 0
    PRINT("final_lambda=",final_lambda,", final_x=",io_util.vector_to_str_n_decimals(final_x,16))
#endif
    #---------------- Check the results. ------------------
    if chart==CH_AKS:
      sigma = kruskal.sigma(x[0],x[1]) # e.g., could have moved from III to II
    if final_lambda>=real_lambda_max:
      BREAK
    x = final_x
    v = final_v
    # Find final Schwarzschild r:
    xs = transform.transform_point(x,spacetime,chart,CH_SCH)
    r = xs[1]
    # Check for incomplete geodesic:
    if r<EPS or opt['dlambda']<EPS:
      r_stuff = runge_kutta.r_stuff(spacetime,chart,x,v,acc,pt,acc_p,pt_p)
      err,r,rdot,rddot,p,lam_left = r_stuff
      final_lambda = final_lambda+lam_left
      err = RK_INCOMPLETE
      BREAK
    opt['triggers'] = user_triggers
    lambda0 = final_lambda
    opt['lambda0'] = final_lambda
  #---------------- Return. ----------------
  if err==RK_INCOMPLETE:
    info['message'] = 'incomplete geodesic'
  return final_helper(err,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,user_chart)

def choose_step_size(r,p,tol,lam_left,x,v,chart):
    """
    Choose step size.
    Returns [n,dlambda,terminate]
    """
    if r<=1.0:
      return choose_step_size_interior(r,p,tol,lam_left)
    else:
      return choose_step_size_exterior(r,tol,x,v,chart)

def choose_step_size_exterior(r,tol,x,v,chart):
    """
    Choose step size.
    Returns [n,dlambda,terminate]
    """
    # Try to estimate an inverse affine-parameter scale for the motion. This is independent of the choice
    # of affine parameter, but in the case of an affine parameter equal to the proper time, this is
    # roughly the inverse time scale for motion by about one schwarzschild radius.
    vs = transform.transform_vector(v,x,SP_SCH,chart,CH_SCH) # find v vector in Sch coords
    scale = abs(vs[0]/r**1.5)+abs(vs[1]/r)+abs(vs[2])+abs(vs[3])+abs(vs[4])
    # The following parameters have to be tuned for optimal performance and only work for
    # a particular order of RK.
    n = 100
    k = 1.0
    order = 4.0 # order of RK
    dlambda = k*tol**(1.0/order)/scale
    return [n,dlambda,FALSE]

def choose_step_size_interior(r,p,tol,lam_left):
    # Use heuristics to pick a step size:
    # fixme -- sanity checks on lam_left
    # fixme: don't hardcode parameters here
    alpha = 0.5
    k = 0.25
    #r_min = 1.0e-16
    r_min = tol**(1.0/p)
    if r_min<1.0e-8:
      r_min=1.0e-8
    # time to quit?
    if r<r_min:
      return [0,0.0,TRUE]
    # set step size
    dlambda = k*tol**0.25*lam_left**alpha
    n=100 # Try to do enough steps with fixed step size to avoid excessive overhead.
    safety = 0.3 # margin of safety so that we never hit singularity
    if n*dlambda>safety*lam_left:
      # ...fixme -- can this be improved?
      n=safety*lam_left/dlambda
      if n<1:
        n=1
        dlambda = safety*lam_left
    return [n,dlambda,FALSE]

def final_helper(err,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,desired_chart):
  x = transform.transform_point(final_x,spacetime,chart,desired_chart)
  v = transform.transform_vector(final_v,final_x,spacetime,chart,desired_chart)
  a = transform.transform_vector(final_a,final_x,spacetime,chart,desired_chart)
  return [err,x,v,a,final_lambda,info,sigma]

def chart_and_triggers(r,triggers,sigma,future_oriented):
  # Pick coordinates that are well adapted to the region we're currently in.
  # Returns optimal chart and creates a list of triggers.
  # fixme -- doesn't work correctly in region III or if not future oriented
  if r>1.1:
    optimal_chart = CH_SCH
    APPEND_TO_ARRAY(triggers,([-1.0,1,1.05,0.3])) # trigger on r<1.05, nearing horizon
  else:
    if r>=1.0:
      optimal_chart = CH_AKS
      # Trigger on a change of sign in either a or b, which means a change of region.
      APPEND_TO_ARRAY(triggers,([ 1.0,0, 0.0,0.3])) # a becoming positive
      APPEND_TO_ARRAY(triggers,([-1.0,0, 0.0,0.3])) # a becoming negative
      APPEND_TO_ARRAY(triggers,([ 1.0,1, 0.0,0.3])) # b becoming positive
      APPEND_TO_ARRAY(triggers,([-1.0,1, 0.0,0.3])) # b becoming negative
    else:
      approaching_singularity = ((sigma>0.0 and future_oriented) or (sigma<0.0 and not future_oriented))
      if r>0.5:
        optimal_chart = CH_SCH
        if approaching_singularity:
          # approaching the singularity
          APPEND_TO_ARRAY(triggers,([-1.0,1,0.35,0.3])) # trigger on r<0.05, nearing singularity
        else:
          # receding from the singularity
          APPEND_TO_ARRAY(triggers,([ 1.0,1,0.95,0.3]))
      else:
        optimal_chart = CH_KEP
        if not approaching_singularity:
          APPEND_TO_ARRAY(triggers,([1.0,1, 0.2,0.3]))
        # It's not useful to try to make a trigger that prevents or detects hitting the singularity. Triggers
        # are too crude for that purpose, don't work reliably because the coordinate velocities diverge.
  return optimal_chart

  
