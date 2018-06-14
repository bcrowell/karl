#include "util.h"
#include "math.h"
#include "io_util.h"
#include "spacetimes.h"
#include "runge_kutta.h"

import schwarzschild,angular

def geodesic_simple(spacetime,chart,x0,v0,opt):
  """
  Calculate a geodesic using geodesic equation and 4th-order Runge-Kutta.

  spacetime = label for the spacetime we're doing (see spacetimes.h for labels)
  chart = label for the coordinate chart we're using
  x = starting point, given as an array of coordinates
  v = components of starting tangent vector (need not be normalized, can be null or spacelike)
  The following options are in the hash opt[].
    lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
                   we hit a singularity)
    dlambda = step size
    ndebug = 0, or, if nonzero, determines how often to print debugging output; e.g., if ndebug=100
    lambda0 = initial affine parameter, defaults to 0
    norm_final = adjust the final x and v to lie on and tangent to the unit sphere in i-j-k space;
                 default=True
  returns
    [err,final_x,final_v,final_lambda,info]
  where
    err = 0 if normal, or bitwise or of codes such as RK_ERR, RK_INCOMPLETE, defined in runge_kutta.h
    final_x,final_v,final_lambda = final values from integration
    info = hash with keys below
      message = error message
  """
  x=copy.deepcopy(x0); v=copy.deepcopy(v0)
  lambda_max=opt["lambda_max"]; dlambda=opt["dlambda"]; ndebug=opt["ndebug"]
  lambda0=0.0
  if hasattr(opt,"lambda0"): lambda0=opt["lambda0"]
  norm_final = True
  if hasattr(opt,"norm_final"): norm_final=opt["norm_final"]
  ok = False
  n = math.ceil(lambda_max/dlambda)
  if ndebug==0:
    steps_between_debugging=n*2 # debugging will never happen
  else:
    steps_between_debugging=ndebug
  debug_count = steps_between_debugging+1 # trigger it on the first iteration
  lam = 0.0
  ok,ndim,christoffel_function = chart_info(spacetime,chart)
  if not ok: return [RK_ERR,x,v,0.0,mess(["unrecognized spacetime or chart: ",spacetime," ",chart])]
  ndim2 = ndim*2 # Reduce 2nd-order ODE to ndim2 coupled 1st-order ODEs.
  if len(x)!=ndim or len(v)!=ndim:
    return [RK_ERR,x,v,0.0,mess(["x or v has wrong length"])]
  order = 4 # 4th order Runge-Kutta
  for iter in range(0,n):
    est = [[0 for i in range(ndim2)] for step in range(order)]
            # =k in the notation of most authors
            # Four estimates of the changes in the independent variables for 4th-order Runge-Kutta.
    debug_count=debug_helper(debug_count,ndebug,steps_between_debugging,iter,lam,x,v)
    y0 = [0 for i in range(ndim2)]
    for i in range(0,ndim): y0[i]=copy.deepcopy(x[i])
    for i in range(0,ndim): y0[i+ndim]=copy.deepcopy(v[i])
    for step in range(0,order):
      if step==0: y=copy.deepcopy(y0)
      if step==1:
        for i in range(0,ndim2):
          y[i] = y0[i]+0.5*est[0][i]
      if step==2:
        for i in range(0,ndim2):
          y[i] = y0[i]+0.5*est[1][i]
      if step==3:
        for i in range(0,ndim2):
          y[i] = y0[i]+est[2][i]
      ch = christoffel_function(y)
      for i in range(0,ndim2): est[step][i]=0.0
      for i in range(0, ndim):
        a = 0.0 # is essentially the acceleration
        for j in range(0, ndim):
          for k in range(0, ndim):
            a -= ch[j][k][i]*y[ndim+j]*y[ndim+k]
        est[step][ndim+i] = a*dlambda
      for i in range(0, ndim):
        est[step][i] = y[ndim+i]*dlambda
    lam= lam+dlambda
    tot_est = [0 for i in range(ndim2)]
    for i in range(0,ndim2):
      tot_est[i] = (est[0][i]+2.0*est[1][i]+2.0*est[2][i]+est[3][i])/6.0
    for i in range(0, ndim):
      v[i] += tot_est[ndim+i]
    for i in range(0, ndim):
      x[i] += tot_est[i]
  debug_helper(debug_count,ndebug,steps_between_debugging,n,lam,x,v)
  if norm_final:
    x = angular.renormalize(x)
    v = angular.make_tangent(x,v)
  return [0,x,v,lam,{}]

def mess(stuff):
  return {'message':strcat(stuff)}

def chart_info(spacetime,chart):
  recognized = False
  if (spacetime|chart)==(SP_SCH|CH_SCH): return [True,5,schwarzschild.christoffel]
  return [False,None,None]

def debug_helper(debug_count,ndebug,steps_between_debugging,iter,lam,x,v):
  """
  Prints debugging info and returns the updated debug_count.
  """
  do_debug = False
  if ndebug!=0 and (debug_count>=steps_between_debugging):
    debug_count = 0
    do_debug = True
  if do_debug:
    print("i=",iter," lam=",("%4.2e" % lam),
                      " x=",io_util.vector_to_str_n_decimals(x,1),
                      " v=",io_util.vector_to_str_n_decimals(v,1))
  return debug_count+1

  
