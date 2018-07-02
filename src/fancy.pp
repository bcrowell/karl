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
  A specialized wrapper for trajectory_simple(), optimized for good numerical precision in the
  Schwarzschild spacetime. Because we may be switching in and out of Schwarzschild or Keplerian
  coordinates, the motion is assumed to be confined to regions I and II.
  Inputs and outputs are the same as for runge_kutta.trajectory_simple().
  """
  if HAS_KEY(opt,'triggers'):
    user_triggers = CLONE_ARRAY_OF_FLOATS2DIM(opt['triggers'])
  else:
    user_triggers = []
  triggers = CLONE_ARRAY_OF_FLOATS2DIM(user_triggers)
  # Find out our initial Schwarzschild r:
  xs = transform.transform_point(x0,spacetime,chart,CH_SCH)
  r = xs[1]
  # Pick coordinates that are well adapted to the region we're currently in.
  if r>1.1:
    optimal_chart = CH_SCH
    APPEND_TO_ARRAY(triggers,([-1.0,1,1.05,0.3])) # trigger on r<1.05, nearing horizon
  else:
    if r>0.9:
      optimal_chart = CH_AKS
      optimal_chart = CH_SCH
      APPEND_TO_ARRAY(triggers,([ 1.0,1, 0.0,0.3])) # trigger on b>0, entering region II
      APPEND_TO_ARRAY(triggers,([-1.0,1, 0.0,0.3])) # trigger on b<0, leaving horizon for region I
    else:
      optimal_chart = CH_KEP
      eps_u = 0.5*r**1.5
      eps_u_min = EPS**1.5 # corresponding to r=EPS
      if eps_u<eps_u_min:
        eps_u=eps_u_min
      APPEND_TO_ARRAY(triggers,([-1.0,1, eps_u,0.1]))
      # ... Trigger on hitting the singularity. Small alpha because it's really bad if we overshoot.
      APPEND_TO_ARRAY(triggers,([1.0,1, 0.926,0.3]))
      # ... Trigger on u>0.926, which is r>0.95. It's not possible for u to trigger this for a causal world-line,
      #     so this would be relevant only for a spacelike geodesic.
  x = transform.transform_point(x0,spacetime,chart,optimal_chart)
  v = transform.transform_vector(v0,x0,spacetime,chart,optimal_chart)
  opt['triggers'] = triggers
  if not HAS_KEY(opt,'debug_chart'):
    opt['debug_chart'] = TRUE
  result = runge_kutta.trajectory_simple(spacetime,optimal_chart,x,v,opt)
  [err,final_x,final_v,final_a,final_lambda,info] = result
  # Find final Schwarzschild r:
  xs = transform.transform_point(x0,spacetime,optimal_chart,CH_SCH)
  r = xs[1]
  # Check for incomplete geodesic:
  if r<10*EPS or opt['dlambda']<10*EPS:
    result[0] = RK_INCOMPLETE
    return result
  opt['triggers'] = user_triggers
  opt['dlambda'] = opt['dlambda']*0.5
  opt['lambda0'] = final_lambda
  return trajectory_schwarzschild(spacetime,optimal_chart,final_x,final_v,opt)
