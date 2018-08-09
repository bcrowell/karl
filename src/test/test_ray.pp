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
  # The following is a smoke test. What tends to happen is that for points very close to the photon sphere,
  # where |a|>>|b| or |b|>>|a|, application of the Christoffel symbols gives a point beyond the
  # singularity, resulting in a crash.
  r = 0.9
  alpha_max = 2.3759564949418355
  d = 1.0e-10
  alpha = alpha_max-d
  beta = alpha_to_beta(r,alpha)

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

r = 0.5 # test inside horizon
test_obs_and_alpha(r)
r = 5.0
test_obs_and_alpha(r)
r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
test_obs_and_alpha(r)

test_alpha_max_schwarzschild()

test_riazuelo_deflection()
test_deflection_continuity()
test_deflection_naughty_cases()
# ... setting time_is_irrelevant in ray.pp fixes this, but breaks other tests
#     maybe it's doing it too often, and errors are accumulating
