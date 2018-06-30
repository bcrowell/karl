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

from io_util import fl

import schwarzschild,kruskal,angular,transform,runge_kutta

def trajectory_schwarzschild(spacetime,chart,x0,v0,opt):
  """
  A specialized version of trajectory_simple(), optimized for good numerical precision in the
  Schwarzschild spacetime. Inputs and outputs are the same as for runge_kutta.trajectory_simple().
  """
  if HAS_KEY(opt,'triggers'):
    user_triggers = CLONE_ARRAY_OF_FLOATS2DIM(opt['triggers'])
  else:
    user_triggers = []
  triggers = CLONE_ARRAY_OF_FLOATS2DIM(user_triggers)
  r_min = x0[1]/100 # assumes Sch. coords; should be /n; should be 1/n if x0[1]>1
  APPEND_TO_ARRAY(triggers,([-1.0,1,r_min,0.3])) # assumes Sch. coords
  opt['triggers'] = triggers
  result = runge_kutta.trajectory_simple(spacetime,chart,x0,v0,opt)
  [err,final_x,final_v,final_a,final_lambda,info] = result
  if final_x[1]<EPS or opt['dlambda']<EPS:
    result[0] = RK_INCOMPLETE
    return result
  opt['triggers'] = user_triggers
  opt['dlambda'] = opt['dlambda']*0.1
  opt['lambda0'] = final_lambda
  return trajectory_schwarzschild(spacetime,chart,final_x,final_v,opt)
