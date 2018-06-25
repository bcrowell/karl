#include "language.h"
#include "util.h"
#include "math.h"
#include "io_util.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#if "LANG" eq "python"
import ctypes,numpy
import c_libs
from c_libs import c_double_p
#...don't import karl_c_lib this way, because you get a copy of it in whatever state it was in?
#    https://stackoverflow.com/a/142601/1142217 -- comment by Paul Whipp
#endif


from io_util import fl

import schwarzschild,kruskal,angular

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
               then we print debugging information at every 100th step
    lambda0 = initial affine parameter, defaults to 0
    norm_final = adjust the final x and v to lie on and tangent to the unit sphere in i-j-k space;
                 default=TRUE
    do_limit_change = boolean, should we do a sanity check by limiting changes in coordinates per step,
            die if limit is violated?; default=FALSE
    limit_change = this is approximately the maximum fractional change in r or 1/10 the maximum change in i,j,k,
            expressed in units of 1/n; default: 1
  returns
    [err,final_x,final_v,final_lambda,info]
  where
    err = 0 if normal, or bitwise or of codes such as RK_ERR, RK_INCOMPLETE, defined in runge_kutta.h
    final_x,final_v,final_lambda = final values from integration
    info = hash with keys below
      message = error message
  """
  x=CLONE_ARRAY_OF_FLOATS(x0)
  v=CLONE_ARRAY_OF_FLOATS(v0)
  lambda_max=opt["lambda_max"]
  dlambda=opt["dlambda"]
  ndebug=opt["ndebug"]
  lambda0=0.0
  if HASATTR(opt,"lambda0"):
    lambda0=opt["lambda0"]
  norm_final = TRUE
  if HASATTR(opt,"norm_final"):
    norm_final=opt["norm_final"]
  n = math.ceil((lambda_max-lambda0)/dlambda)
  do_limit_change = FALSE
  if HASATTR(opt,"do_limit_change"):
    do_limit_change=opt["do_limit_change"]
  if do_limit_change:
    limit_change = 1.0/n
    if HASATTR(opt,"limit_change"):
      do_limit_change=float(opt["do_limit_change"])/n
  ok = FALSE
  if ndebug==0:
    steps_between_debugging=n*2 # debugging will never happen
  else:
    steps_between_debugging=ndebug
  debug_count = steps_between_debugging+1 # trigger it on the first iteration
  lam = lambda0
  ok,ndim,christoffel_function = chart_info(spacetime,chart)
#if "LANG" eq "js"
  use_c = false; __NO_TRANSLATION__
#endif
#if "LANG" eq "python"
  use_c = ((spacetime|chart)==(SP_SCH|CH_SCH))
  use_c = False # doesn't work yet
#endif
  if not ok:
    return [RK_ERR,x,v,0.0,mess(["unrecognized spacetime or chart: ",spacetime," ",chart])]
  ndim2 = ndim*2 # Reduce 2nd-order ODE to ndim2 coupled 1st-order ODEs.
  if LEN(x)!=ndim or LEN(v)!=ndim:
    return [RK_ERR,x,v,0.0,mess(["x or v has wrong length"])]
  order = 4 # 4th order Runge-Kutta
#if "LANG" eq "python"
  # Allocate some arrays that are later repeatedly reused.
  acc = numpy.zeros((ndim,))
  acc = acc.astype(numpy.float64)
  pt = numpy.zeros((ndim2,))
  pt = pt.astype(numpy.float64)
#endif
  for iter in range(0,n):
    est = [[0 for i in range(ndim2)] for step in range(order)] #js est=karl.array2d(ndim2,order);
    #         =k in the notation of most authors
    #         Four estimates of the changes in the independent variables for 4th-order Runge-Kutta.
    debug_count=debug_helper(debug_count,ndebug,steps_between_debugging,iter,lam,x,v)
    y0 = EMPTY1DIM(ndim2)
    for i in range(0,ndim):
      y0[i]=CLONE_FLOAT(x[i])
    for i in range(0,ndim):
      y0[i+ndim]=CLONE_FLOAT(v[i])
    for step in range(0,order):
      if step==0:
        y=CLONE_ARRAY_OF_FLOATS(y0)
      if step==1:
        for i in range(0,ndim2):
          y[i] = y0[i]+0.5*est[0][i]
      if step==2:
        for i in range(0,ndim2):
          y[i] = y0[i]+0.5*est[1][i]
      if step==3:
        for i in range(0,ndim2):
          y[i] = y0[i]+est[2][i]
      for i in range(0,ndim2): est[step][i]=0.0
      #use_c = False
      if use_c:
        # use faster C implementation:
#if "LANG" eq "python"
        for i in range(0, ndim2):
          p[i]=y[i]
        c_libs.karl_c_lib.apply_christoffel(spacetime,chart,
                pt.ctypes.data_as(c_double_p),
                acc.ctypes.data_as(c_double_p),
                ctypes.c_double(dlambda))
        for i in range(0, ndim):
          est[step][ndim+i] = acc[i] # C routine takes care of multiplying a by dlambda
#endif
      else:
        ch = christoffel_function(y)
        for i in range(0, ndim):
          a = 0.0 # is essentially the acceleration
          for j in range(0, ndim):
            for k in range(0, ndim):
              a -= ch[j][k][i]*y[ndim+j]*y[ndim+k]
          est[step][ndim+i] = a*dlambda
      for i in range(0, ndim):
        est[step][i] = y[ndim+i]*dlambda
    lam= lam+dlambda
    tot_est = EMPTY1DIM(ndim2)
    for i in range(0,ndim2):
      tot_est[i] = (est[0][i]+2.0*est[1][i]+2.0*est[2][i]+est[3][i])/6.0
    for i in range(0, ndim):
      v[i] += tot_est[ndim+i]
    for i in range(0, ndim):
      if do_limit_change:
        check_limit_change(spacetime,chart,x,tot_est,limit_change)
      x[i] += tot_est[i]
  debug_helper(debug_count,ndebug,steps_between_debugging,n,lam,x,v)
  if norm_final:
    x = angular.renormalize(x)
    v = angular.make_tangent(x,v)
  return [0,x,v,lam,{}]

def check_limit_change(spacetime,chart,x,dx,limit_change):
  """
  Sanity check to flag sudden large changes in coordinates.
  """
  ok = TRUE
  if (spacetime|chart)==(SP_SCH|CH_SCH):
    rel_dr=abs(dx[1])/x[1]
  if (spacetime|chart)==(SP_SCH|CH_AKS):
    rel_dr=(abs(dx[0])+abs(dx[1]))/(1+abs(x[0]-x[1]))
    # ... quick and dirty estimate using r=a-b+1, not really appropriate for small distances
  if rel_dr>limit_change:
    THROW(io_util.strcat(['r changed by too much , rel_dr=',rel_dr, \
                                                     ', x=',io_util.vector_to_str(x), \
                                                     ', dx=',io_util.vector_to_str(dx)]))
  for i in range(2,5):
    if abs(dx[i])>10.0*limit_change:
      THROW('angular coord. changed by too much')

def mess(stuff):
  return {'message':io_util.strcat(stuff)}

def chart_info(spacetime,chart):
  recognized = FALSE
  if (spacetime|chart)==(SP_SCH|CH_SCH):
    return [TRUE,5,schwarzschild.christoffel]
  if (spacetime|chart)==(SP_SCH|CH_AKS):
    return [TRUE,5,kruskal.christoffel]
  return [FALSE,None,None]

def debug_helper(debug_count,ndebug,steps_between_debugging,iter,lam,x,v):
  """
  Prints debugging info and returns the updated debug_count.
  """
  do_debug = FALSE
  if ndebug!=0 and (debug_count>=steps_between_debugging):
    debug_count = 0
    do_debug = TRUE
  if do_debug:
    PRINT("i=",iter," lam=",io_util.fl(lam), \
                      " x=",io_util.vector_to_str_n_decimals(x,1), \
                      " v=",io_util.vector_to_str_n_decimals(v,1))
  return debug_count+1

  
