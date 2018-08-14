#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"
#include "spacetimes.h"

import test,vector,ray,star_properties,transform
import scipy.integrate as integrate

verbosity = 1

def test_alpha_max_schwarzschild():
  assert_equal_eps(ray.alpha_max_schwarzschild(15.0),3.02,0.02)
  # ... estimating asymptote from Riazuelo fig. 1
  assert_equal_eps(ray.alpha_max_schwarzschild(1.0001),MATH_PI-acos(23/31),1.0e-4)
  # ... Riazuelo eq. 86 (note that his 82.1 deg. should be 84.2 deg.)
  assert_equal_eps(ray.alpha_max_schwarzschild(1.0001),ray.alpha_max_schwarzschild(0.9999),1.0e-4)
  # ... should be a continuous function as we cross the horizon
  assert_equal_eps(ray.alpha_max_schwarzschild(1.50001),ray.alpha_max_schwarzschild(1.49999),1.0e-5)
  # ... logic changes at r=3/2, but result should be a continuous function of r
  assert_equal_eps(ray.alpha_max_schwarzschild(1.5),acos(-sqrt(2/3)),1.0e-7)
  # ... I don't actually know that this is correct, but it interpolates in a sensible way between the
  #     known-good results at r=1 and r=infty, and it's what I seem to get from the code. The conjectured
  #     value smells right, since the photon sphere involves expressions with sqrt(3) in them.
  #     The value in radians is 2.52611294491941.
  assert_equal_eps(ray.alpha_max_schwarzschild(1.0e-16),MATH_PI/2,1.0e-7)
  # ... Riazuelo, p. 15, eq. 87

def test_obs_and_alpha(r):
  test_schwarzschild_standard_observer(r) 
  test_le_to_alpha_schwarzschild(r)
  test_alpha_to_le_schwarzschild(r)

def test_alpha_to_le_schwarzschild(r):
  # Test the trivial special cases of alpha=0 and pi:
  le,in_n_out = ray.alpha_to_le_schwarzschild(0,r)
  assert_equal_eps(le,0.0,10*EPS)
  le,in_n_out = ray.alpha_to_le_schwarzschild(MATH_PI,r)
  assert_equal_eps(le,0.0,1.0e-14)
  # Test for correct round-trip behavior through the function and its inverse:
  for i in range(10):
    alpha = 0.1+(2.9/9)*i # values of alpha from 0.1 to 3.0
    test_alpha_to_le_schwarzschild_round_trip(r,alpha)

def test_alpha_to_le_schwarzschild_round_trip(r,alpha):
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x_obs,v_obs,rho,j = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  le,in_n_out = ray.alpha_to_le_schwarzschild(alpha,r)
  alpha2,vv = ray.le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  #print("r=",r,", alpha=",alpha,", le=",le,", in_n_out=",in_n_out,", alpha2=",alpha2)
  assert_equal_eps(alpha2,alpha,1.0e-12)

def test_deflection_naughty_cases():
  # If these are going to fail, they basically tend to fail by crashing or hanging up forever, i.e., these
  # are essentially smoke tests. But if L/E is nonzero, we also compare with the exact result from
  # reducing the problem to quadrature. (Could also do this for L/E=0, but would have to use a different
  # equation.)
  # What tends to happen is that for points very close to the photon sphere,
  # where |a|>>|b| or |b|>>|a|, application of the Christoffel symbols gives a point beyond the
  # singularity, resulting in a crash. The thing that seems to fix all of these for r>~0.1 is that in
  # do_ray_schwarzschild2(), I check whether norm is nonzero, and if it is, I retry with smaller
  # lambda. Before I figured out that adaptive algorithm, I was also seeing cases where there
  # was no crash into the singularity, but inaccurate results came out (e.g., there would be
  # an obvious glitch in the aberration table, where beta(alpha) would not be monotonic).
  # Some of the crashing cases below were found by playing around with values of alpha close
  # to such a glitchy value. I also have logic now in make_aberration_table where if there is
  # a glitch in the aberration table, it tries to detect and patch it, and it prints out a warning.
  #---
  eps = 1.0e-6 # default tolerance compared to quadrature (relative error)
  test_number = 1
  #---
  # Small initial r value, seems to require step size propto r^2 rather than r.
  r = 0.1
  alpha = 1.0
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,eps)
  #---
  # Chokes at the very end, when it gets to large r. Norm creeps up bigger and bigger, finally
  # just barely crosses the threshold of 1.0e-6 so that it triggers an error.
  r = 0.39681846091970896
  alpha = 2.085642127611168 # close to alpha_max=2.1526252824191947
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,eps)
  #---
  # A trajectory parallel to the horizon. If we use kruskal_to_time_zero() too enthusiastically, this
  # trajectory gets sucked toward (a,b)=0 and ends up with nonzero norm.
  r = 0.9
  alpha = 0.0
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,eps)
  #---
  # If dlambda isn't small enough, crashes with bad norm near horizon.
  r = 0.6080841337780387
  alpha = 2.200625558361262 # close to alpha_max=2.2646182578916436
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,eps)
  #---
  # Requires small step size for values of r up to about 1.01. Otherwise it
  # crashes for these inputs, and also gives inaccurate results for nearby values of alpha.
  r = 0.608084133778
  alpha = 0.8387
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,eps)
  #---
  # A demanding case where a past-oriented geodesic spends a long time in the photon sphere before
  # emerging. Doesn't have as small an error as the others.
  r = 0.9
  d = 1.0e-10
  alpha = ray.alpha_max_schwarzschild(r)-d
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,0.005)
  #---
  # Passing the following test seems to require a very small step size in do_ray_schwarzschild2().
  # See logic involving dlambda_safety.
  # With larger step sizes, what seems to happen is that as we pass out through the horizon, the
  # norm of the velocity vector becomes unacceptably large (-0.02) due to numerical errors.
  r = 0.5863
  alpha = 0.7085
  test_number = test_deflection_naughty_one(verbosity,r,alpha,test_number,eps)

def test_deflection_naughty_one(verbosity,r,alpha,test_number,eps):
  le = ray.alpha_to_le_schwarzschild(alpha,r)[0]
  if verbosity>=2:
    print("testing naughty case ",test_number,", r=",r,", alpha=",alpha,", |L/E|=",le)
  beta = alpha_to_beta(r,alpha)
  if le!=0.0:
    p = le**-2 # notation used by Gibbons, arxiv.org/abs/1110.6508
    q = integrate.quad(lambda r: 1/sqrt(p*r**4-r**2+r), r, numpy.inf)[0]
    assert_rel_equal_eps(beta,q,eps)
  if verbosity>=2:
    print("  beta=",beta)
    if le!=0.0:
      print("  result from quadrature=",q," rel. err.=",(beta-q)/q)
  return test_number+1

def test_deflection_continuity():
  # Currently, my algorithms have qualitative differences for r in (0,1), (1,1.5), and (1.5,infty).
  # Make sure that deflection angles appear to be continuous across these boundaries.
  # In the following, the alphas are close to the maximums, to make it a severe test.
  # --
  beta1 = alpha_to_beta(1.001,2.3)
  beta2 = alpha_to_beta(0.999,2.3)
  assert_equal_eps(beta1,beta2,0.01)
  # --
  beta1 = alpha_to_beta(1.50001,2.5)
  beta2 = alpha_to_beta(1.49999,2.5)
  assert_equal_eps(beta1,beta2,0.001)

def test_riazuelo_deflection():
  r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
#if "LANG" eq "python"
  # Interpolating from his graph, he has lines crossing at alpha=beta=2.48+-0.1.
  alpha = 2.48
  beta = alpha_to_beta(r,alpha)
  assert_equal_eps(alpha,2.48,0.1)
  assert_equal_eps(beta,alpha,0.1)

  # Check sign of aberration for outward angles, should have alpha>beta because SR dominates:
  alpha = 0.1
  beta = alpha_to_beta(r,alpha)
  assert_boolean(alpha>beta,"should have alpha>beta at outward angles, because SR dominates")  

  # Check sign of aberration for grazing angles, should have alpha<beta because deflection dominates:
  alpha = 2.9
  beta = alpha_to_beta(r,alpha)
  assert_boolean(alpha<beta,"should have alpha<beta at grazing angles, because deflection dominates")  
#endif

def alpha_to_beta(r,alpha):
  tol = 1.0e-6
  count_winding(0.0,[],[],0,0,{})
  beta,done,final_v,region = ray.do_ray_schwarzschild(r,tol,count_winding,alpha)
  return beta

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


def test_le_to_alpha_schwarzschild(r):
  le_over_max = 0.1
  alpha,v = test_le_to_alpha_schwarzschild_one(r,le_over_max,1) # glancing inward

  alpha,v = test_le_to_alpha_schwarzschild_one(r,0.0,0) # directly outward
  assert_equal_eps(alpha,0.0,1.0e-6)

  alpha,v = test_le_to_alpha_schwarzschild_one(r,0.0,1) # directly inward
  assert_equal_eps(alpha,MATH_PI,1.0e-6)


def test_le_to_alpha_schwarzschild_one(r,le_over_max,in_n_out):
  # le_over_max = |L/E|/|L/E|_max
  #    For r<1 there really is no max value of L/E, but we arbitrarily set the max to r/sqrt|1-1/r|.
  #    See https://physics.stackexchange.com/q/418157/4552 .
  # in_n_out = 0 for outward radial motion
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x_obs,v_obs,rho,j = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1-1/r
  le_max = r/sqrt(abs(aa)) # abs is so we still get a test at r<1, see above
  le = le_over_max*le_max
  alpha,v = ray.le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  # v has norm 0:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,v),0.0,10*EPS)
  # v is in future light cone:
  is_future,a_plus_b = transform.sch_is_in_future_light_cone(x_obs,v)
  if a_plus_b<=0.0:
    PRINT("r=",r,", a_plus_b=",a_plus_b,", v=",v)
  assert_boolean(a_plus_b>0.0,"v of photon is not future-oriented")

  return [alpha,v]

def test_schwarzschild_standard_observer(r):
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}

  x_obs,v_obs,rho,j = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1.0-1.0/r

  # observer's energy is 1:
  assert_equal_eps(aa*v_obs[0],1.0,10*EPS) 

  # rho and j both have norm -1:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,rho),-1.0,10*EPS)
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,j),-1.0,10*EPS)

  # v_obs has norm 1:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,v_obs),1.0,10*EPS)

  # v_obs, rho, and j are all orthogonal to one another:
  assert_equal_eps(vector.inner_product(spacetime,chart,pars,x_obs,v_obs,rho),0.0,10*EPS)
  assert_equal_eps(vector.inner_product(spacetime,chart,pars,x_obs,v_obs,j),0.0,10*EPS)
  assert_equal_eps(vector.inner_product(spacetime,chart,pars,x_obs,rho,j),0.0,10*EPS)

  # rho has a positive r component:
  if r>1.0:
    assert_boolean(rho[1]>0.0,"observer's radial vector has a negative r component in Schwarzschild coordinates")

  # j has a positive j component:
  assert_boolean(j[3]>0.0,"observer's azimuthal vector j has a negative j component")

  # rho has vanishing angular components, observer is supposed to be radially infalling:
  assert_equal_eps(rho[2],0.0,10*EPS)
  assert_equal_eps(rho[3],0.0,10*EPS)
  assert_equal_eps(rho[4],0.0,10*EPS)

  # v_obs is oriented in the future-timelike direction:
  is_future,a_plus_b = transform.sch_is_in_future_light_cone(x_obs,v_obs)
  assert_boolean(is_future,"v_obs is not future-oriented")

test_deflection_naughty_cases()

r = 0.5 # test inside horizon
test_obs_and_alpha(r)
r = 5.0
test_obs_and_alpha(r)
r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
test_obs_and_alpha(r)

test_alpha_max_schwarzschild()

test_riazuelo_deflection()
test_deflection_continuity()
