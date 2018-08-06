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
       star_properties

#--------------------------------------------------------------------------------------------------

#if "LANG" eq "python"
# ... fixme: js doesn't work with the way I did the nested sub and closure for count_winding()

def make_aberration_tables(r,tol,verbosity):
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
  last_deflection = 0.0
  s0 = 4
  # ...scaling factor for number of angular steps
  #    Making this a big value, like 10, makes fake stars take a long time, because the sky is subdivided
  #    very finely. Making it a very small value, like 2, may make interpolation of Doppler shifts too crude,
  #    and also causes poor behavior for r just inside the photon sphere.
  s1 = 2 # for moderate deflections, reduce step size by this factor
  s2 = 10 # ... additional factor for large deflections
  n_angles = s0*s1*s2 # we actually skip most of these angles
  for in_n_out in range(2): # 0=outward, 1=inward
    for i in range(n_angles):
      frac_skip = s1*s2 # normally we fill in most of the angles by interpolation
      if last_deflection>0.5:
        frac_skip = frac_skip//s1
      if last_deflection>2.0:
        frac_skip = frac_skip//s2
      skip_this = not (i%frac_skip==0)
      z = (float(i)/float(n_angles)) # note that this should *not* be n_angles-1
      if in_n_out==0:
        le = z*max_le
      else:
        if is_schwarzschild:
          le_ph = 0.5*3**1.5 # L/E of the photon sphere, unstable circular orbits for photons in Schwarzschild
          le = math_util.linear_interp(0.0,1.0,max_le,le_ph,z)
        else:
          # more generic logic, more likely to be appropriate if not Sch.
          le = (1.0-z)*max_le
      alpha,v_observation = le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
      x = x_obs # initial position
      v = CLONE_ARRAY_OF_FLOATS(v_observation)
      # ... Initial velocity of ray, which is simulated going backward in time, as if emitted rather than absorbed.
      #     Will get rescaled later, see comments below, but clone it to make sure it can't get munged.
      if skip_this:
        beta,done = [NONE,FALSE] # fill in later by interpolation
      else:
        beta,done,v_emission = do_ray(spacetime,chart,pars,x,v,r,tol,count_winding,alpha)
      table.append([r,alpha,beta])
      got_result = not (IS_NONE(beta) or IS_NAN(beta))
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
      if got_result:
        q = v_emission[0]
        v_observation = vector.scalar_mult(v_observation,1.0/q)
        # This vector was set up for an outgoing ray, but we're flipping this around, treating it as an incoming
        # ray that was observed here. Therefore we need to flip components 1...4:
        v_observation = vector.scalar_mult(v_observation,-1.0)
        v_observation[0] = -v_observation[0] # flip 0 component back to what it was
        v_table.append([alpha,]+v_observation)
        # ... The + here is concatenation of the lists. This won't work in js.
      if got_result:
        if verbosity>=2:
          pass
          #PRINT("r=",r,", L/E=",le,", alpha=",alpha*180.0/MATH_PI," deg, beta=",beta*180.0/MATH_PI," deg.")
        if verbosity>=3:
          if alpha==0.0:
            err_approx = 0.0
          else:
            approx = -0.5*alpha/(sqrt(r)-1)
            err_approx = (approx-abs(alpha-beta))/alpha
            print("alpha=",alpha,", beta-alpha=",beta-alpha,", approx=",approx)
            #PRINT("r=",r,", alpha=",alpha*180.0/MATH_PI," deg, beta=",beta*180.0/MATH_PI,\
            #    " deg, rel err in approx=",err_approx)
        last_deflection = abs(alpha-beta)
        if abs(alpha-beta)>5.0*MATH_PI: # Riazuelo says 5pi is enough to get all visual effects.
          PRINT("Deflection=",abs(alpha-beta)*180.0/MATH_PI," deg. is >5pi, done.")
          done = TRUE
      if done:
        break
    if done:
      break
  fill_in_aberration_table_by_interpolation(table,r)
  table2 = []
  for i in range(len(table)):
    x = CLONE_ARRAY_OF_FLOATS(table[i])
    # x = [r,alpha,beta]
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
    x.append(beta-alpha)
    x.append(f)
    table2.append(x)
  return [table2,v_table]

#endif

def do_ray(spacetime,chart,pars,x,v,r,tol,count_winding,alpha):
  n = 100
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  ri = r
  lambda_max = 10.0 # fixme, sort of random
  info = {}
  WHILE(TRUE)
    opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,'tol':tol,\
          'user_function':count_winding,'no_enter_horizon':TRUE}
    err,final_x,final_v,final_a,final_lambda,info,sigma  = \
            fancy.trajectory_schwarzschild(spacetime,chart,pars,x,v,opt)
    if err!=0:
      if err==RK_INCOMPLETE or err==RK_TRIGGER:
        PRINT("Geodesic at alpha=",alpha*180.0/MATH_PI," deg. is incomplete or entered horizon, done.")
        return [NAN,TRUE,NONE]
      THROW("error: "+str(err))
    rf = final_x[1]
    if IS_NAN(rf):
      THROW("rf is NaN")
    if rf>1.0e9 and rf>100.0*r:
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
  return [beta,FALSE,final_v]

def fill_in_aberration_table_by_interpolation(table,r):
  #-------------------
  # Use interpolation to fill in missing values.
  # Find i and j that both have real data.
  for i in range(len(table)-2):
    if IS_NONE(table[i][2]):
      continue
    j=i+1
    WHILE(j<=len(table)-1 and IS_NONE(table[j][2]))
      j=j+1
    END_WHILE
    if j==i+1 or j>len(table)-1:
      continue
    alpha1 = table[i][1]
    alpha2 = table[j][1]
    beta1 = table[i][2]
    beta2 = table[j][2]
    # Find error compared to small-angle approx:
    err1 = alpha1/(sqrt(r)-1)-(beta1-alpha1)
    err2 = alpha2/(sqrt(r)-1)-(beta2-alpha2)
    k=i+1
    WHILE(IS_NONE(table[k][2]) and k<j)
      alpha = table[k][1]
      err = math_util.linear_interp(alpha1,alpha2,err1,err2,alpha) # interpolate to find est. of what error would have been
      beta = alpha+alpha/(sqrt(r)-1)-err
      table[k][2] = beta
      k=k+1
    END_WHILE
  # Remove any uncomputed items from the end of the table:
  WHILE(len(table)>0 and (IS_NONE(table[len(table)-1][2]) or IS_NAN(table[len(table)-1][2])))
    POP_FROM_ARRAY(table)
  END_WHILE
  if len(table)==0:
    THROW("entire array removed")

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
  if not transform.sch_is_in_future_light_cone(x_obs,v):
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


