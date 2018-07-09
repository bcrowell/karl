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

import schwarzschild,kruskal,keplerian,angular,transform

def trajectory_simple(spacetime,chart,x0,v0,opt):
  """
  Calculate a trajectory using geodesic equation plus external force term, with 4th-order Runge-Kutta.

  spacetime = label for the spacetime we're doing (see spacetimes.h for labels)
  chart = label for the coordinate chart we're using
  x = starting point, given as an array of coordinates
  v = components of starting tangent vector
    The velocity need not be normalized, can be null or spacelike. If it is normalized, then the
    affine parameter represents proper time.
  The following options are in the hash opt[].
    lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
                   we hit a singularity)
    dlambda = step size
    ndebug = 0, or, if nonzero, determines how often to print debugging output; e.g., if ndebug=100
               then we print debugging information at every 100th step
    debug_function = a function to be called when printing debugging output
    lambda0 = initial affine parameter, defaults to 0
    norm_final = adjust the final x and v to lie on and tangent to the unit sphere in i-j-k space;
                 default=TRUE
    triggers = array of 4-element arrays, each describing a trigger (see below)
    force_acts = boolean, do we have an external force?
    force_function = function that calculates the proper acceleration vector d^2x/dlambda^2,
                             given (lambda,x,v) as inputs; its output will be used and immediately discarded,
                             so the function does not need to clone it before returning it
    force_chart = chart that the function wants for its inputs and outputs
  triggers
    These allow the integration to be halted when it appears that in the next iteration,
    a certain coordinate or velocity would cross a certain threshold.
    [0] = sense, +1 or -1 for triggering in a rising or falling direction
    [1] = index of coordinate (0-4) or velocity (5-9) on which to trigger
    [2] = threshold value
    [3] = fudge factor alpha, which should be less than 1 if something bad happens at threshold or if it would
                   be bad not to get the trigger
  returns
    [err,final_x,final_v,final_a,final_lambda,info]
  where
    err = 0 if normal, or bitwise or of codes such as RK_ERR, RK_INCOMPLETE, defined in runge_kutta.h
    final_x,final_v,final_a,final_lambda = final values of position, velocity, acceleration, and affine param
    info = hash with keys below
      message = error message
  """
  x=CLONE_ARRAY_OF_FLOATS(x0)
  v=CLONE_ARRAY_OF_FLOATS(v0)
  #-- process input options
  lambda_max,dlambda,ndebug,debug_function,lambda0,norm_final,\
        n_triggers,trigger_s,trigger_on,trigger_threshold,trigger_alpha,\
        force_acts,force_function,force_chart =\
        runge_kutta_get_options_helper(opt)
  #-- initial setup
  n,steps_between_debugging,debug_count,lam,ok,ndim,christoffel_function = \
           runge_kutta_init_helper(lambda_max,lambda0,dlambda,ndebug,spacetime,chart)
  if not ok:
    return [RK_ERR,x,v,0.0,mess(["unrecognized spacetime or chart: ",spacetime," ",chart])]
  use_c = c_available(spacetime,chart)
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
  acc_p = acc.ctypes.data_as(c_double_p)
  pt_p = pt.ctypes.data_as(c_double_p)
#else
  acc = EMPTY1DIM(ndim)
#endif
  y0 = EMPTY1DIM(ndim2)
  for iter in range(n):
    dlambda = (lambda_max-lam)/(n-iter) # small readjustment so we land on the right final lambda
    est = [[0 for i in range(ndim2)] for step in range(order)] #js est=karl.array2d(ndim2,order);
    #         =k in the notation of most authors
    #         Four estimates of the changes in the independent variables for 4th-order Runge-Kutta.
    debug_count=debug_helper(debug_count,ndebug,steps_between_debugging,iter,lam,dlambda,\
                                        x,v,debug_function,spacetime|chart)
    for i in range(ndim):
      y0[i]=x[i]
    for i in range(ndim):
      y0[i+ndim]=v[i]
    y0 = CLONE_ARRAY_OF_FLOATS(y0)
    # ...Disentangle it from x and v so that in python, changing x or v can't change it. This is actually
    #    not necessary, because x and v don't change until the next iteration, when y0 is built again,
    #    but I find it too hard to reason about the code without this.
    for step in range(order):
      if step==0:
        y=CLONE_ARRAY_OF_FLOATS(y0)
      if step==1:
        for i in range(ndim2):
          y[i] = y0[i]+0.5*est[0][i]
      if step==2:
        for i in range(ndim2):
          y[i] = y0[i]+0.5*est[1][i]
      if step==3:
        for i in range(ndim2):
          y[i] = y0[i]+est[2][i]
      for i in range(ndim2):
        est[step][i]=0.0
      if use_c:
        # use faster C implementation:
#if "LANG" eq "python"
        for i in range(ndim2):
          pt[i]=y[i]
        c_libs.karl_c_lib.apply_christoffel(spacetime,chart,pt_p,acc_p,ctypes.c_double(dlambda))
#endif
      else:
        apply_christoffel(christoffel_function,y,acc,dlambda,ndim)
      if force_acts:
        handle_force(acc,lam,x,v,force_function,force_chart,ndim,spacetime,chart,dlambda)
      for i in range(ndim):
        est[step][ndim+i] = acc[i]
        est[step][i] = y[ndim+i]*dlambda
    if n_triggers>0 and \
            trigger_helper(x,v,acc,dlambda,n_triggers,trigger_s,trigger_on,trigger_threshold,trigger_alpha,ndim):
      return runge_kutta_final_helper(debug_count,ndebug,steps_between_debugging,iter,lam,dlambda,x,v,acc,\
                     norm_final,debug_function,chart)
    #-- Update everything:
    lam= lam+dlambda
    tot_est = EMPTY1DIM(ndim2)
    for i in range(ndim2):
      tot_est[i] = (est[0][i]+2.0*est[1][i]+2.0*est[2][i]+est[3][i])/6.0
    for i in range(ndim):
      v[i] += tot_est[ndim+i]
    for i in range(ndim):
      x[i] += tot_est[i]
  return runge_kutta_final_helper(debug_count,ndebug,steps_between_debugging,n,lam,dlambda,x,v,acc,norm_final,\
              debug_function,spacetime|chart)

def c_available(spacetime,chart):
#if "LANG" eq "js"
  return false; __NO_TRANSLATION__
#endif
#if "LANG" eq "python"
  return (spacetime==SP_SCH and (chart==CH_SCH or chart==CH_AKS or chart==CH_KEP))
#endif

def r_stuff(spacetime,chart,x,v,acc,pt,acc_p,pt_p):
  """
  Returns [r,r',r'',p,lam_left], where the primes represent derivatives with respect to affine parameter,
  p is an estimate of the exponent in r ~ lambda^p, and lam_left is an estimate of the distance
  left before the singularity in terms of the affine parameter.

  The arrays acc_p and pt_p are pointers to arrays that have already been allocated, or
  None if this is the js implementation.
  """
  ndim = 5
  ndim2 = 10
  x2 = transform.transform_point(x,spacetime,chart,CH_SCH)
  v2 = transform.transform_vector(v,x,spacetime,chart,CH_SCH)
  for i in range(ndim):
    pt[i] = x2[i]
    pt[i+5] = v2[i]
  if c_available(SP_SCH,CH_SCH):
    # use faster C implementation:
#if "LANG" eq "python"
    c_libs.karl_c_lib.apply_christoffel(SP_SCH,CH_SCH,pt_p,acc_p,ctypes.c_double(1.0))
#endif
  else:
    ok,ndim,christoffel_function,name = transform.chart_info(SP_SCH|CH_SCH)
    apply_christoffel(christoffel_function,pt,acc,1.0,ndim)
  r = x2[1]
  rdot = v2[1]
  rddot = acc[1]
  p = 1.0/(1-r*rddot/(rdot*rdot))
  if p<0.4:
    p=0.4
  if p>1.0:
    p=1.0
  lam_left = -p*r/rdot # estimate of when we'd hit the singularity
  return [r,rdot,rddot,p,lam_left]

def handle_force(a,lam,x,v,force_function,force_chart,ndim,spacetime,chart,dlambda):
  # The API says that force_function does not need to clone its output vector before returning it, so
  # we need to make sure to discard it here and never do anything with it later.
  # The vector a that we're modifying already has a factor of dlambda in it, so we multiply by dlambda here
  # as well.
  x2 = transform.transform_point(x,spacetime,chart,force_chart)
  v2 = transform.transform_vector(v,x,spacetime,chart,force_chart)
  proper_accel2 = force_function(lam,x2,v2)
  proper_accel = transform.transform_vector(proper_accel2,x2,spacetime,force_chart,chart)
  for i in range(ndim):
    a[i] = a[i]+proper_accel[i]*dlambda

def runge_kutta_get_options_helper(opt):
  lambda_max  =runge_kutta_get_par_helper(opt,"lambda_max",NONE)
  dlambda     =runge_kutta_get_par_helper(opt,"dlambda",NONE)
  ndebug      =runge_kutta_get_par_helper(opt,"ndebug",0)
  debug_function =runge_kutta_get_par_helper(opt,"debug_function",default_debug_function)
  lambda0     =runge_kutta_get_par_helper(opt,"lambda0",0.0)
  norm_final  =runge_kutta_get_par_helper(opt,"norm_final",TRUE)
  n_triggers = 0
  trigger_s,trigger_on,trigger_threshold,trigger_alpha = [[],[],[],[]]
  if HAS_KEY(opt,"triggers"):
    n_triggers = runge_kutta_get_trigger_options_helper(opt,trigger_s,trigger_on,trigger_threshold,trigger_alpha)
  force_acts     =runge_kutta_get_par_helper(opt,"force_acts",FALSE)
  force_function =runge_kutta_get_par_helper(opt,"force_function",0)
  force_chart     =runge_kutta_get_par_helper(opt,"force_chart",0)
  return [lambda_max,dlambda,ndebug,debug_function,lambda0,norm_final, \
                   n_triggers,trigger_s,trigger_on,trigger_threshold,trigger_alpha, \
                   force_acts,force_function,force_chart]

def runge_kutta_init_helper(lambda_max,lambda0,dlambda,ndebug,spacetime,chart):
  n = CEIL((lambda_max-lambda0)/dlambda) # dlambda will be adjusted slightly in order to deal with the rounding
  if ndebug==0:
    steps_between_debugging=n*2 # debugging will never happen
  else:
    steps_between_debugging=ndebug
  debug_count = steps_between_debugging+1 # trigger it on the first iteration
  lam = lambda0
  ok,ndim,christoffel_function,name = transform.chart_info(spacetime|chart)
  return [n,steps_between_debugging,debug_count,lam,ok,ndim,christoffel_function]

def trigger_helper(x,v,acc,dlambda,n_triggers,trigger_s,trigger_on,trigger_threshold,trigger_alpha,ndim):
  """
  Check triggers:
  We can trigger in the rising (s=+1) or falling (s=-1) direction. The coordinate or velocity
  we're triggering on differs from the trigger value by dx, and it's currently changing at
  a rate x_dot. Depending on the signs of s, dx, and x_dot, we have 8 cases. The logic below
  handles all the cases properly. The input variable acc is actually the acceleration times dlambda.
  """
  for i in range(n_triggers):
    s = trigger_s[i] # sense of the trigger (see above)
    m = trigger_on[i] # index of coordinate or velocity on which to trigger
    thr = trigger_threshold[i] # threshold value
    alpha = trigger_alpha[i] # fudge factor, can typically be 1; see docs
    if m<ndim:
      # triggering on a coordinate
      dx = thr-x[m]
      x_dot = v[m]
    else:
      # triggering on a velocity
      dx = thr-v[m-ndim]
      x_dot = acc[m-ndim]/dlambda # left over from step==3, good enough for an estimate
    #PRINT("s=",s,", dx=",dx,", x_dot=",x_dot,", dlambda=",dlambda,"lhs=",x_dot*dlambda*s,", rhs=",alpha*dx*s)
    if s*dx>0 and x_dot*dlambda*s>alpha*dx*s: # Note that we can't cancel the s, which may be negative.
      # We extrapolate that if we were to complete this iteration, we would cross the threshold.
      return TRUE
  return FALSE

def runge_kutta_final_helper(debug_count,ndebug,steps_between_debugging,n,lam,dlambda,x,v,acc,norm_final,debug_function,spacetime_and_chart):
  debug_helper(debug_count,ndebug,steps_between_debugging,n,lam,dlambda,x,v,debug_function,spacetime_and_chart)
  # ... always do a printout for the final iteratation
  if norm_final:
    x = angular.renormalize(x)
    v = angular.make_tangent(x,v)
  return [0,x,v,acc,lam,{}]

def runge_kutta_get_par_helper(opt,name,default_val):
  val = default_val
  if HAS_KEY(opt,name):
    val = opt[name]
  if IS_NONE(val):
    THROW('required option '+name+' not supplied')
  return val

def runge_kutta_get_trigger_options_helper(opt,trigger_s,trigger_on,trigger_threshold,trigger_alpha):
  n_triggers = LEN(opt["triggers"])
  for i in range(n_triggers):
    trigger = opt["triggers"][i]
    APPEND_TO_ARRAY(trigger_s,trigger[0])
    APPEND_TO_ARRAY(trigger_on,trigger[1])
    APPEND_TO_ARRAY(trigger_threshold,trigger[2])
    APPEND_TO_ARRAY(trigger_alpha,trigger[3])
  return n_triggers

def apply_christoffel(christoffel_function,y,acc,dlambda,ndim):
  # A similar routine, written in C for speed, is in apply_christoffel.c. This version exists
  # only so it can be translated into javascript.
  ch = christoffel_function(y)
  for i in range(ndim):
    a = 0.0 # is essentially the acceleration
    for j in range(ndim):
      for k in range(ndim):
        a -= ch[j][k][i]*y[ndim+j]*y[ndim+k]
    acc[i] = a*dlambda

def mess(stuff):
  return {'message':io_util.strcat(stuff)}

def debug_helper(debug_count,ndebug,steps_between_debugging,iter,lam,dlambda,\
                                x,v,debug_function,spacetime_and_chart):
  """
  Prints debugging info and returns the updated debug_count.
  """
  do_debug = FALSE
  if ndebug!=0 and (debug_count>=steps_between_debugging):
    debug_count = 0
    do_debug = TRUE
  if do_debug:
    ok,ndim,christoffel_function,name = transform.chart_info(spacetime_and_chart) # just to get name of chart
    debug_function(iter,lam,dlambda,x,v,name)
  return debug_count+1

def default_debug_function(iter,lam,dlambda,x,v,name):
  PRINT("i=",iter," lam=",io_util.fl(lam), \
                      " x=",io_util.vector_to_str_n_decimals(x,1), \
                      " v=",io_util.vector_to_str_n_decimals(v,1))
  
