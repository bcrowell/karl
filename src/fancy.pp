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

def trajectory_schwarzschild(spacetime,chart,x0,v0,opt,sigma):
  """
  A specialized wrapper for trajectory_simple(), optimized for speed and numerical precision in the
  Schwarzschild spacetime.
  Inputs and outputs are the same as for runge_kutta.trajectory_simple(), plus the additional input sigma.

  Values of spacetime and chart are defined in spacetimes.h.
  x0 and v0 = initial position and velocity, expressed in the chosen chart
  opt = same hash of options defined in runge_kutta
  sigma = see docs, distinguishes regions I and II from III and IV, if this ambiguity isn't resolvable
             in the chosen coordinate chart
  returns [err,x,v,a,final_lambda,info,sigma]
  """
  if HAS_KEY(opt,'triggers'):
    user_triggers = CLONE_ARRAY_OF_FLOATS2DIM(opt['triggers'])
  else:
    user_triggers = []
  n_unproductive_iterations = 0
  lambda0 = 0.0
  if HAS_KEY(opt,'lambda0'):
    lambda0 = opt['lambda0']
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
          eps_u = 0.5*r**1.5 # trigger on u becoming half its present value
          eps_u_min = EPS**1.5 # corresponding to r=EPS
          if eps_u<eps_u_min:
            eps_u=eps_u_min
          APPEND_TO_ARRAY(triggers,([-1.0,1, eps_u,0.1]))
          # ... Trigger on hitting the singularity. Small alpha because it's really bad if we overshoot.
          APPEND_TO_ARRAY(triggers,([1.0,1, 0.2,0.3]))
          # ... only relevant for a spacelike world-line
    if chart!=optimal_chart:
      x2 = transform.transform_point(x,spacetime,chart,optimal_chart)
      v2 = transform.transform_vector(v,x,spacetime,chart,optimal_chart)
      x = x2
      v = v2
    chart = optimal_chart
    opt['triggers'] = triggers
    result = runge_kutta.trajectory_simple(spacetime,chart,x,v,opt)
    err,final_x,final_v,final_a,final_lambda,info = result
    if chart==CH_AKS:
      sigma = kruskal.sigma(x[0],x[1]) # e.g., could have moved from III to II
    if final_lambda>=opt['lambda_max']:
      return final_helper(err,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,user_chart)
    # Count consecutive iterations in which we didn't make any progress in terms of lambda:
    if final_lambda==lambda0:
      n_unproductive_iterations = n_unproductive_iterations+1
    else:
      n_unproductive_iterations = 0
    if n_unproductive_iterations>20:
      return final_helper(RK_ERR,final_x,final_v,final_a,final_lambda,\
                          {'message':'too many unproductive iterations'},\
                          sigma,spacetime,chart,user_chart)
    x = final_x
    v = final_v
    # Find final Schwarzschild r:
    xs = transform.transform_point(x,spacetime,chart,CH_SCH)
    r = xs[1]
    # Check for incomplete geodesic:
    if r<10*EPS or opt['dlambda']<10*EPS:
      result[0] = RK_INCOMPLETE
      return final_helper(RK_INCOMPLETE,final_x,final_v,final_a,final_lambda,\
                          {'message':'incomplete geodesic'},\
                          sigma,spacetime,chart,user_chart)
    opt['triggers'] = user_triggers
    lambda0 = final_lambda
    opt['lambda0'] = final_lambda
    if chart==CH_KEP and r<1.0: # popped out of RK because we were about to hit singularity on the next step
      opt['dlambda'] = opt['dlambda']*0.5

def final_helper(err,final_x,final_v,final_a,final_lambda,info,sigma,spacetime,chart,desired_chart):
  x = transform.transform_point(final_x,spacetime,chart,desired_chart)
  v = transform.transform_vector(final_v,final_x,spacetime,chart,desired_chart)
  a = transform.transform_vector(final_a,final_x,spacetime,chart,desired_chart)
  return [err,x,v,a,final_lambda,info,sigma]

