# This module contains general-relativistic ray-tracing code that isn't
# related to human vision or rendering, doesn't need to read or write
# image files or databases, and doesn't use the python ephem library.
# Everything in this module is supposed to translate properly to js, but
# some functions currently don't, and are commented out with conditionals.

#include "language.h"
#include "util.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild,euclidean,celestial,math_util,\
       star_properties,kruskal
import scipy.integrate as integrate

#--------------------------------------------------------------------------------------------------

#if "LANG" eq "python"
# ... fixme: js doesn't work with the way I did the nested sub and closure for count_winding()

def make_aberration_tables(r,tol,verbosity,max_deflection):
  """
  Determine a table of optical aberration angles for an observer in the Schwarzschild spacetime.
  Each line of the table is in the format [r,alpha,beta,beta-alpha,f].
  Here r is observer's coordinate
  in units of Schwarzschild radius, alpha is angle as seen by observer if they're in a standard
  state of motion (free fall from rest at infinity), beta is angle on the celestial sphere,
  and f is the factor by which intensities are amplified by lensing. 
  In addition, we return a separate table [alpha,v0,...v4] of the 
  components of the velocity vector v^kappa of the light ray at observation,
  in 5-dimensional Schwarzschild coordinates. 
  This is used in order to find the Doppler shifts later.
  The alpha values in this table are a subset of the
  ones in the aberration table, because it seems safer to interpolate the Doppler shifts at the end,
  rather than interpolating the vector components here.
  The normalization of the velocities is such that the velocity vector at emission is (1,-1,0,0,0).
  The observer's velocity is not calculated here because it can easily be found later, and
  would be the same for all rays: -- u_kappa=(1,-A*sqrt(1-A)), where A=1-1/r. Note that this
  is the covariant form of u.
  """
  is_schwarzschild = TRUE
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  x_obs,v_obs,rho,j = schwarzschild_standard_observer(r,spacetime,chart,pars)
  if r>1.0:
    max_le = r/sqrt(aa) # maximum possible L/E for a photon at this r (not the max. for visibility)
  else:
    max_le = 10.0 # fixme
  done = FALSE
  def count_winding(lam,x,v,spacetime,chart,pars):
    if lam==0.0:
      count_winding.winding = 0
      count_winding.last_angle = 0.0
    else:
      angle = atan2(x[3],x[2]) # returns an angle from -pi to pi
      if angle<0.0:
        angle = angle+2.0*MATH_PI
      if count_winding.last_angle>5.0 and angle<1.0:
        count_winding.winding = count_winding.winding+1
      if count_winding.last_angle<1.0 and angle>5.0:
        count_winding.winding = count_winding.winding-1
      count_winding.last_angle = angle
    return count_winding.winding
  count_winding(0.0,[],[],0,0,{})
  table = []
  v_table = []
  n_angles = 100
  alpha_max = alpha_max_schwarzschild(r)
  print("  Calculating table of aberration angles, alpha_max=",alpha_max,", ",n_angles," angles")
  for i in range(n_angles):
    ii = (n_angles-1)-i # ii decreases, so that z decreases and alpha increases
    z = (float(ii+1)/float(n_angles)) # z varies from 1/n_angles to 1
    alpha = alpha_max*(1-z**4) # points are closely spaced near the vertical asymptote
#if 0
    print("--------------------------------------------------------------------------------------")
    print("r=",r,", alpha=",alpha)
#endif
    le,in_n_out = alpha_to_le_schwarzschild(alpha,r)
    alpha2,v_observation = le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
    x = x_obs # initial position
    v = CLONE_ARRAY_OF_FLOATS(v_observation)
    # ... Initial velocity of ray, which is simulated going backward in time, as if emitted rather than absorbed.
    #     Will get rescaled later, see comments below, but clone it to make sure it can't get munged.
    d = alpha_max-alpha
    beta,done,v_emission,region = do_ray_schwarzschild(r,tol,count_winding,alpha,0)
    if done:
      BREAK # can get incomplete geodesic at alpha<alpha_max due to numerical precision
    if abs(beta-alpha)>max_deflection:
      PRINT("  Quitting at deflection ",abs(beta-alpha),", greater than limit of ",max_deflection)
      BREAK
    table.append([r,alpha,beta,beta-alpha])
    smooth = beta+log(d)
    #print("r=",r,", alpha=",alpha,", beta=",beta,", region=",region[2])
    #---
    # Calculate the velocity of the ray at observation. This would actually be pretty trivial, since
    # v is the "initial" velocity for solving the diffeq, but is actually the final
    # velocity of the ray, at observation. However, we do need v_emission for normalization.
    # Because the Schwarzschild coordinates are spherical coordinates, v_emission is guaranteed
    # to have the form v^kappa=(q,-q,0,0,0). We need to retrofit v=v_observation to have a normalization
    # that would have resulted from emission with v^kappa=(1,-1,0,0,0), which corresponds to a different  
    # choice of affine parameter. In the Schwarzschild case, v_t is conserved, so the result for
    # v_observation should always be v^t=1/(1-1/r), and we would not need the results of ray-tracing
    # to tell us this, but I want this code to be more general, so I don't use that assumption.
    q = v_emission[0]
    v_observation = vector.scalar_mult(v_observation,1.0/q)
    # This vector was set up for an outgoing ray, but we're flipping this around, treating it as an incoming
    # ray that was observed here. Therefore we need to flip components 1...4:
    v_observation = vector.scalar_mult(v_observation,-1.0)
    v_observation[0] = -v_observation[0] # flip 0 component back to what it was
    v_table.append([alpha,]+v_observation)
    # ... The + here is concatenation of the lists. This won't work in js.
  if len(table)==0:
    THROW('no data points in aberration table')
  table2 = []
  # Check for glitches:
  for i in range(len(table)-2):
    r,alpha1,beta1,diff = table[i]
    r,alpha2,beta2,diff = table[i+1]
    r,alpha3,beta3,diff = table[i+2]
    if (beta2-beta1)*(beta3-beta2)<0.0:
      # slope changes sign, shouldn't happen because beta(alpha) is monotonic
      PRINT("warning: apparent glitch repaired in aberration table, r=",r,", alpha=",alpha2,", beta=",beta2)
      table[i+1][2] = 0.5*(beta1+beta3)
  for i in range(len(table)):
    x = CLONE_ARRAY_OF_FLOATS(table[i])
    # x = [r,alpha,beta,beta-alpha]
    if i==0:
      ii=0
      jj=1
    else:
      ii=i-1
      jj=i
    alpha = 0.5*(table[ii][1]+table[jj][1]) # average them because we don't want alpha/beta=0/0 at i=0
    beta  = 0.5*(table[ii][2]+table[jj][2])
    dalpha = table[jj][1]-table[ii][1]
    dbeta  = table[jj][2]-table[ii][2]
    f = abs(sin(alpha)*dalpha/(sin(beta)*dbeta))
    if verbosity>=2:
      pass
      #print("alpha=",alpha,", beta=",beta,", dalpha=",dalpha,", dbeta=",dbeta,", f=",f)
    # =dOmega(obs)/dOmega(infinity)=amplification, by Liouville's thm (?)
    # is abs() right?; beta can be >pi
    x.append(f)
    table2.append(x)
  return [table2,v_table]

#endif

def alpha_to_beta_by_quadrature(r,alpha):
  """
  In the Schwarzschild spacetime, determine the deflection of a ray arriving to the standard observer
  from null infinity, using numerical integration. Doesn't check whether alpha is legal.
  Cf. alpha_to_beta().
  """
  le = alpha_to_le_schwarzschild(alpha,r)[0]
  if le==0.0:
    return 0.0 # is the correct result, if we haven't been supplied with an illegal alpha=pi
  p = le**-2 # notation used by Gibbons, arxiv.org/abs/1110.6508
  result = integrate.quad(lambda r: 1/sqrt(p*r**4-r**2+r), r, numpy.inf)
  q = result[0]
  if result[1]>1.0e-4:
    THROW('high error estimate')
  return q

def alpha_to_beta(r,alpha):
  """
  In the Schwarzschild spacetime, determine the deflection of a ray arriving to the standard observer
  from null infinity, using numerical integration.
  Cf. alpha_to_beta_by_quadrature().
  """
  def count_winding(lam,x,v,spacetime,chart,pars):
    if lam==0.0:
      count_winding.winding = 0
      count_winding.last_angle = 0.0
    else:
      angle = atan2(x[3],x[2]) # returns an angle from -pi to pi
      if angle<0.0:
        angle = angle+2.0*MATH_PI
      if count_winding.last_angle>5.0 and angle<1.0:
        count_winding.winding = count_winding.winding+1
      if count_winding.last_angle<1.0 and angle>5.0:
        count_winding.winding = count_winding.winding-1
      count_winding.last_angle = angle
    return count_winding.winding
  tol = 1.0e-6
  count_winding(0.0,[],[],0,0,{})
  beta,done,final_v,region = do_ray_schwarzschild(r,tol,count_winding,alpha,0)
  return beta

def do_ray_schwarzschild(r,tol,count_winding,alpha,force_algorithm):
  if r<0.3:
    algorithm=3
  else:
    if r<1.0:
      algorithm=2
    else:
      algorithm=1
  if force_algorithm!=0:
    algorithm=force_algorithm
  if algorithm==1:
    return do_ray_schwarzschild1(r,tol,count_winding,alpha)
  if algorithm==2:
    return do_ray_schwarzschild2(r,tol,count_winding,alpha)
  if algorithm==3:
    return do_ray_schwarzschild3(r,tol,count_winding,alpha)

def do_ray_schwarzschild1(r,tol,count_winding,alpha):
  """
  Returns [beta,if_incomplete,final_v,region].
  This algorithm tends to fail for past-oriented geodesics that exit the event horizon. It oscillates
  back and forth across the horizon, I think because of either a bug or a lack of numerical stability.
  The alternative do_ray_schwarzschild2() is meant to work for those cases.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x,v_obs,rho,j = schwarzschild_standard_observer(r,spacetime,chart,pars)
  le,in_n_out = alpha_to_le_schwarzschild(alpha,r)
  alpha2,v = le_to_alpha_schwarzschild(r,le,in_n_out,x,v_obs,rho,spacetime,chart,pars)
  # don't need alpha, just need v of photon
  n = 100
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  ri = r
  r_obs = r
  lambda_max = 10.0 # fixme, sort of random
  info = {}
  lam = 0.0
  if r>1.0:
    no_enter_horizon = TRUE
  else:
    no_enter_horizon = FALSE
  WHILE(TRUE)
#if 0
    print("x=",io_util.vector_to_str(x),", v=",io_util.vector_to_str(v))
#endif
#if 0
    print("before fancy, r=",x[1])
#endif
    opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,'tol':tol,\
          'user_function':count_winding,'no_enter_horizon':no_enter_horizon}
    err,final_x,final_v,final_a,final_lambda,info,sigma  = \
            fancy.trajectory_schwarzschild(spacetime,chart,pars,x,v,opt)
#if 0
    print("after fancy, r=",final_x[1])
#endif
    lam = lam+final_lambda
    if err!=0:
      if err==RK_INCOMPLETE or err==RK_TRIGGER:
        PRINT("Geodesic at alpha=",alpha," radians is incomplete or entered horizon, err=",err,", done.")
        return [NAN,TRUE,NONE,[2,0,"interior"]]
      THROW("error: "+str(err))
    rf = final_x[1]
    if IS_NAN(rf):
      THROW("rf is NaN")
    if rf>1.0e9*r_obs and rf>100.0*r:
      break
    if rf>10.0 and final_v[1]>0 and rf>ri: # quick and dirty test, far away and getting farther
      f = (rf-ri)/ri # fractional change in r from the iteration we just completed
      if f<0.1:
        lambda_max = lambda_max*(0.1/f) # try to get at least a 10% change in r with each iteration
    x = final_x
    v = final_v
    ri = rf
  END_WHILE
  w = info['user_data'] # winding number
  beta = atan2(final_x[3],final_x[2]) # returns an angle from -pi to pi
  if beta<0.0:
    beta = beta+2.0*MATH_PI
  beta = beta+w*2.0*MATH_PI
  return [beta,FALSE,final_v,[1,0,"exterior"]]

def do_ray_schwarzschild2(r,tol,count_winding,alpha):
  """
  Returns [beta,if_incomplete,final_v,region].
  Algorithm that only uses KS coordinates, assumes we're starting inside
  horizon and heading outward. Calls do_ray_schwarzschild2_one_try(), and normally
  that works and we're done. But if the "bad" flag is set, we start over from scratch
  with a smaller basic step size and see if that helps.
  The parameter tol is normally 1.0e-6;
  this is not really calibrated against anything absolute, but has an effect on both the
  step size and the maximum norm that we allow.
  """
  n = 100
  max_norm = tol
  basic_dlambda = 0.01*(tol/1.0e-6)**0.25
  max_retries_from_scratch = 2
  n_retries_from_scratch = 0
  WHILE(n_retries_from_scratch<max_retries_from_scratch)
#if 0
    if n_retries_from_scratch>0:
      PRINT("retrying from scratch, basic_dlambda=",basic_dlambda)
#endif
    result = do_ray_schwarzschild2_one_try(r,basic_dlambda,count_winding,alpha,n,max_norm)
    retry_code = result.pop()
    if retry_code==0:
      return result
    if retry_code==2:
      THROW('failed with norm too big after first block')
    n_retries_from_scratch = n_retries_from_scratch+1
    basic_dlambda = basic_dlambda*0.1
  END_WHILE

def do_ray_schwarzschild2_one_try(r,basic_dlambda,count_winding,alpha,n,max_norm):
  """
  Returns [beta,if_incomplete,final_v,region,retry_code].
  Algorithm that only uses KS coordinates, assumes we're starting inside
  horizon and heading outward. Do Runge-Kutta in blocks of n iterations.
  If the norm gets bigger than max_norm, indicating
  an inaccurate result, we first try redoing the most recent block with
  exponentially decreasing step sizes. If that doesn't work, we return
  with retry_code=1. A normal return is retry_code=0, and retry_code=2 means
  that we failed on the first block, so retrying won't help.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x,v_obs,rho,j = schwarzschild_standard_observer(r,spacetime,chart,pars)
  le,in_n_out = alpha_to_le_schwarzschild(alpha,r)
  alpha2,v = le_to_alpha_schwarzschild(r,le,in_n_out,x,v_obs,rho,spacetime,chart,pars)
  # don't need alpha, just need v of photon
  #----
  # Convert from Schwarzschild to KS.
  x_ks = transform.transform_point(x,spacetime,chart,pars,CH_AKS)
  v_ks = transform.transform_vector(v,x,spacetime,chart,pars,CH_AKS)
  chart = CH_AKS
  x = x_ks
  v = v_ks
  # Flip a and b components of velocity so it's past-oriented. Don't flip angular part, because
  # that has the effect of flipping the sign of beta.
  v[0] = -v[0]
  v[1] = -v[1]
  #----
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  ri = r
  lambda_max = 0.2 # fixme, sort of random
  dlambda = basic_dlambda
  info = {}
  lam = 0.0
#define DEBUG_DO_RAY_SCHWARZSCHILD2 0
#if DEBUG_DO_RAY_SCHWARZSCHILD2
  print("\n\n\n\n\n\n\n\n\n\n=============================================================\n\n\n\n\n\n\n\n\n\n")
  print("  starting, observer at r=",r,", alpha=",alpha)
#endif
  max_retries = 3
  block = 0
  WHILE(TRUE)
    block = block+1
    dlambda_safety = 1.0
    if ri>1.0:
      dlambda_small_r_safety = 1.0
    else:
      dlambda_small_r_safety = ri*ri
    dlambda_actual = dlambda*dlambda_safety*dlambda_small_r_safety
    opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,\
          'user_function':count_winding,'dlambda':dlambda_actual,'ndebug':0,'time_is_irrelevant':TRUE}
#if DEBUG_DO_RAY_SCHWARZSCHILD2
    print("before simple")
    print("  r=",ri,", lambda=",lam," dlambda=",dlambda," norm=",\
           vector.norm(spacetime,chart,pars,x,v))
    print("  ri=",ri," dlambda_actual=",dlambda_actual)
    print("  x=",io_util.vector_to_str(x),", v=",io_util.vector_to_str(v))
#endif
    n_retry = 0
    WHILE(TRUE)
      err,final_x,final_v,final_a,final_lambda,info  = \
              runge_kutta.trajectory_simple(spacetime,chart,pars,x,v,opt)
      tf,rf,mu = kruskal.aux(final_x[0],final_x[1])
      norm = vector.norm(spacetime,chart,pars,final_x,final_v)
      if abs(norm)<max_norm:
        if n_retry>0 and verbosity>=2:
          PRINT("norm=",norm," is OK after ",n_retry," retries")
        BREAK #--- normal exit from loop ---
      if verbosity>=2 or (DEBUG_DO_RAY_SCHWARZSCHILD2==1):
        PRINT("norm of geodesic is ",norm,", not close to zero, ri=",ri,", rf=",rf,", alpha=",alpha,\
               ",dlambda=",dlambda_actual)
      n_retry = n_retry+1
      if n_retry>max_retries:
        if block>1:
          retry_code = 1
        else:
          # In the situation where it fails on the first block, redoing from scratch won't help.
          retry_code = 2
        return [NAN,FALSE,NONE,kruskal.describe_region(-final_x[0],-final_x[1]),retry_code]
      dlambda_safety = dlambda_safety*0.1
      if verbosity>=2 or (DEBUG_DO_RAY_SCHWARZSCHILD2==1):
        PRINT("*** retry ",n_retry,": redoing with safety=",dlambda_safety)
      dlambda_actual = dlambda*dlambda_safety*dlambda_small_r_safety
      opt['dlambda']=dlambda_actual
    END_WHILE
#if DEBUG_DO_RAY_SCHWARZSCHILD2
    if rf<5.0:
      print("after simple, r=",rf)
#endif
#if DEBUG_DO_RAY_SCHWARZSCHILD2
    if ri<1.0 and rf>1.0:
      print("crossed horizon")
    if ri<1.5 and rf>1.5:
      print("crossed photon sphere")
    if ri<5.0 and rf>5.0:
      print("crossed r=5")
#endif
    lam = lam+final_lambda
    if err!=0:
      if err==RK_INCOMPLETE:
        PRINT("Geodesic at alpha=",alpha," radians is incomplete, done.")
        return [NAN,TRUE,NONE,kruskal.describe_region(-x[0],-x[1]),0]
      THROW("error: "+str(err))
    if IS_NAN(rf):
      THROW("rf is NaN")
    if rf>1.0e6 and rf>100.0*r:
      #print("quitting, rf=",rf,", lambda=",lam,", dlambda=",dlambda)
      break
    if rf>1.7 and rf>ri: # make sure we're well clear of the photon sphere, where rf-ri may be small
      f = (rf-ri)/ri # fractional change in r from the iteration we just completed
      desired_f = 0.2
      if rf>5.0:
        desired_f = 1.0
      if rf>20.0:
        desired_f = 3.0
      if f<desired_f:
        z = (desired_f/f) # try to get at least a certain fractional change in r with each iteration
        lambda_max = lambda_max*z
        dlambda = dlambda*z
    x = final_x
    v = final_v
    ri = rf
  END_WHILE
  w = info['user_data'] # winding number
  beta = atan2(final_x[3],final_x[2]) # returns an angle from -pi to pi
  if beta<0.0:
    beta = beta+2.0*MATH_PI
  beta = beta+w*2.0*MATH_PI
  return [beta,FALSE,final_v,kruskal.describe_region(-x[0],-x[1]),0]
  # ...Minus signs because something about the way I'm treating the time-reversed ray tracing seems
  #    to result in exiting into III when it should be I. Presumably rays that exit into I are
  #    really exiting into III...? fixme -- investigate this

def do_ray_schwarzschild3(r,tol,count_winding,alpha):
  """
  Calculate deflection of a ray by using the fact that null geodesic motion can be reduced to quadrature.
  Returns [beta,if_incomplete,final_v,region].
  The good news is that this is very efficient and reliable, even at small r.
  The bad news is that this type of algorithm will not work for most other spacetimes.
  """
  alpha_max = alpha_max_schwarzschild(r)
  if abs(alpha)>alpha_max:
    return [NAN,TRUE,[],[]]
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  x_obs,v_obs,rho,j = schwarzschild_standard_observer(r,spacetime,chart,pars)
  le,in_n_out = alpha_to_le_schwarzschild(alpha,r)
  beta = alpha_to_beta_by_quadrature(r,alpha)
  alpha2,v = le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  # don't need alpha, just need v of photon
  final_v_t = v[0]*aa # Av_t is conserved in Schwarzschild spacetime, and A=1 at infinity.
  final_v_r = final_v_t # Schwarzschild metric has g_tt=g_rr at infinity.
  final_v = [final_v_t,final_v_r,0.0,0.0,0.0]
  return [beta,FALSE,final_v,[1,0,"exterior"]]

def schwarzschild_standard_observer(r,spacetime,chart,pars):
  """
  The main input is the Schwarzschild radial coordinate r of an observer in the Schwarzschild spacetime.
  The other inputs are described in comments in the calling code.
  Although the inputs include data about the spacetime and chart, this code will not actually work except
  for Schwarzschild cooordinates in the Schwarzschild spacetime.
  The output is the four vectors [x_obs,v_obs,rho,j] describing an observer at the standard state of motion,
  free-falling from rest at infinity. These are the observer's position and velocity;
  the observer's radial vector rho, parallel-transported in from infinity; and the observer's azimuthal
  unit vector (in the direction of increasing j coordinate in the 5-dimensional coordinates).
  """
  aa = 1-1/r
  x_obs = [0.0,r,1.0,0.0,0.0]
  # Calculate the velocity vector of the observer, for the standard state of motion:
  v_obs = [1/aa,-sqrt(1-aa),0.0,0.0,0.0]
  #    ... solution of the equations E=Atdot=1, |v|^2=Atdot^2-(1/A)rdot^2=1
  #        We check in test_ray that it actually is a solution.
  # Calculate the observer's radial vector rho, parallel-transported in from infinity.
  rho = vector.proj(spacetime,chart,pars,x_obs,v_obs,[0.0,1.0,0.0,0.0,0.0])
  # ... take rhat and project out the part orthogonal to the observer's velocity; not yet normalized
  rho = vector.normalize_spacelike(spacetime,chart,pars,x_obs,rho)
  # ... now normalize it to have |rho|^2=-1
  # Gram-Schmidt process to determine j:
  j = vector.proj(spacetime,chart,pars,x_obs,v_obs,[0.0,0.0,0.0,1.0,0.0])
  j = vector.proj_spacelike(spacetime,chart,pars,x_obs,rho,j)
  j = vector.normalize_spacelike(spacetime,chart,pars,x_obs,j)
  # Various checks on these vectors are done in test_ray.
  return [x_obs,v_obs,rho,j]

def alpha_max_schwarzschild(r):
  """
  Maximum alpha at which an observer at r, in the standard state of motion, receives rays.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  le_ph = 0.5*3**1.5 # L/E of the photon sphere, unstable circular orbits for photons in Schwarzschild
  x_obs,v_obs,rho,j = schwarzschild_standard_observer(r,spacetime,chart,pars)
  if r>=1.5:
    in_n_out = 1 # ray peeled off of the photon sphere and came outward to us
  else:
    in_n_out = 0 # limiting rays peel off of the photon sphere and fall inward to us
  if r<1.0:
    in_n_out = 1
    # ...I'm not clear on why this is necessary, but the opposite choice gives discontinuous
    #    and obviously wrong results, while this choice gives a correct result at r=1-epsilon.
  alpha_max,v_observation = le_to_alpha_schwarzschild(r,le_ph,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  return alpha_max

def alpha_to_le_schwarzschild(alpha,r):
  """
  Inputs are Schwarzschild r coordinate and alpha, the angle away from the zenith as measured by the
  standard observer. Return value is [|L/E|,in_n_out]. 
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x_obs,v_obs,rho,j = schwarzschild_standard_observer(r,spacetime,chart,pars)
  v_rho = vector.scalar_mult(rho,cos(alpha))
  v_j   = vector.scalar_mult(j,sin(alpha))
  v_s = vector.add(v_rho,v_j) # purely spacelike part of photon's velocity vector; is normalized by construction
  aa = 1-1/r
  # Loop over four possibilities and pick the one that gives the closest result on the inverse transformation:
  smallest_error = 999.9
  for timelike_sign in range(2):
    # Either of the following gives v.v=0.
    if timelike_sign==0:
      v = vector.add(v_s,v_obs)
    else:
      v = vector.sub(v_s,v_obs)
    le = abs((r*r*v[3])/(aa*v[0]))
    for in_n_out in range(2):
      alpha2,vv = le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
      err = abs(alpha2-alpha)
      if err<smallest_error:
        best = [le,in_n_out]
        smallest_error = err
  return best

def le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars):
  """
  The main inputs are Schwarzschild radius r of the observer, the ratio L/E of the angular momentum of a photon
  to its energy, and in_n_out, which equals 0 if the time-reversed photon is going radially outward, 1 if inward.
  The sign of L/E is ignored.
  The other inputs are described in comments in the calling code.
  Although the inputs include data about the spacetime and chart, this code will not actually work except
  for Schwarzschild cooordinates in the Schwarzschild spacetime.
  Output is [alpha,v], where v is the velocity vector of the photon at observation
  and alpha is the azimuthal angle of the photon as measured by the observer.
  """
  aa = 1-1/r
  #----
  # Find velocity vector of the photon.
  if le==0.0:
    if in_n_out==0:
      v_r = aa
    else:
      # Normally this case is not of much interest, because we don't receive photons along the line from the b.h.
      v_r = -aa
    v = [1.0,v_r,0.0,0.0,0.0]
  else:
    z = r*r/(le*le)-aa
    if z<EPS and z>-10.0*EPS:
      z=EPS # tiny neg z happens due to rounding; make z big enough so it doesn't cause div by zero below
    if z<0.0:
      THROW("z<0, a photon with angular momentum this big can't exist at this r")
      # ... won't happen with current algorithm
    dphi_dr = (1.0/r)/sqrt(z)
    dphi_dt = abs(le*aa)/(r*r)
    dr_dt = abs(dphi_dt/dphi_dr)
    if in_n_out==1:
      dr_dt = -dr_dt
    v = [1.0,dr_dt,0.0,dphi_dt,0.0] # tangent vector
  #----
  is_future,a_plus_b = transform.sch_is_in_future_light_cone(x_obs,v)
  if not (a_plus_b>0):
    v = vector.scalar_mult(v,-1.0)

  alpha = v_photon_and_obs_frame_to_alpha(spacetime,chart,pars,x_obs,v,v_obs,rho)
  return [alpha,v]

def v_photon_and_obs_frame_to_alpha(spacetime,chart,pars,x_obs,v,v_obs,rho):
  v_perp = vector.proj(spacetime,chart,pars,x_obs,v_obs,do_time_refl_schwarzschild(v))
  # ... the part of the photon's velocity orthogonal to the observer's velocity.
  #     The time reflection is necessary empirically in order to get the right direction for
  #     the kinematic aberration, and sort of makes sense because we're doing reverse ray
  #     tracing, but I'm uneasy about the underlying logic and whether this makes sense on the interior.
  v_perp = vector.normalize_spacelike(spacetime,chart,pars,x_obs,v_perp)
  # ... normalized
  zzz = -vector.inner_product(spacetime,chart,pars,x_obs,rho,v_perp)
  # ... Minus sign is because this is +--- signature, so euclidean
  #     dot products have flipped signs.
  alpha = acos(math_util.force_into_range(zzz,-1.0,1.0))
  # ... angle at which the observer says the photon is emitted, see docs; the force_into_range() is
  #     necessary because sometimes we get values for zzz like 1.0000000000000002 due to rounding
  return alpha

def do_time_refl_schwarzschild(v0):
  """
  Reflect the vector spatially in the static frame. The vector is assumed to be in Schwarzschild coordinates.
  Only makes sense for r>1, there is no static frame for r<=1.
  """
  v = CLONE_ARRAY_OF_FLOATS(v0)
  v[0] = -v[0]
  return v


