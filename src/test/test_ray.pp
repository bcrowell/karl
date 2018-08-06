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

def test_riazuelo_deflection():
  r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
#if "LANG" eq "python"
  # Interpolating from his graph, he has lines crossing at alpha=beta=2.48+-0.1.
  # Working backwards, this requires about L/E=11.3, which is 0.728 of max:
  #     calc -e "r=15; aa=1-1/r; lemax=r/sqrtaa; le=11.3/lemax"
  # This is an important test of whether I'm getting SR aberration and GR deflection right;
  # this is a case in which the two effects cancel out. Previously I had been getting this
  # wrong because of confusion about hendling of aberration with time-reversed ray tracing.
  le_over_max = 0.728
  in_n_out = 1
  alpha,beta = test_riazuelo_deflection_one(r,le_over_max,in_n_out)
  assert_equal_eps(alpha,2.48,0.1)
  assert_equal_eps(beta,alpha,0.1)

  # Check sign of aberration for outward angles, should have alpha>beta because SR dominates:
  le_over_max = 0.1
  in_n_out = 0
  alpha,beta = test_riazuelo_deflection_one(r,le_over_max,in_n_out)
  assert_boolean(alpha>beta,"should have alpha>beta at outward angles, because SR dominates")  

  # Check sign of aberration for grazing angles, should have alpha<beta because deflection dominates:
  le_over_max = 0.2 # anything large compared to 1/((3/2)r)=0.04 should avoid absorption; 0.2 gives big defl
  in_n_out = 1
  alpha,beta = test_riazuelo_deflection_one(r,le_over_max,in_n_out)
  assert_boolean(alpha<beta,"should have alpha<beta at grazing angles, because deflection dominates")  
#endif

def test_riazuelo_deflection_one(r,le_over_max,in_n_out):
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x_obs,v_obs,rho = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1-1/r
  le_max = r/sqrt(abs(aa)) # abs is so we still get a test at r<1, see above
  le = le_over_max*le_max
  alpha,v = ray.le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  tol = 1.0e-6
  count_winding(0.0,[],[],0,0,{})
  beta,done,final_v = ray.do_ray(spacetime,chart,pars,x_obs,v,r,tol,count_winding,alpha)
  return [alpha,beta]

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
  x_obs,v_obs,rho = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1-1/r
  le_max = r/sqrt(abs(aa)) # abs is so we still get a test at r<1, see above
  le = le_over_max*le_max
  alpha,v = ray.le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  # v has norm 0:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,v),0.0,10*EPS)
  # v is in future light cone:
  assert_boolean(transform.sch_is_in_future_light_cone(x_obs,v),"v of photon is not future-oriented")

  return [alpha,v]

def test_schwarzschild_standard_observer(r):
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}

  x_obs,v_obs,rho = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1.0-1.0/r

  # observer's energy is 1:
  assert_equal_eps(aa*v_obs[0],1.0,10*EPS) 

  # rho has norm -1:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,rho),-1.0,10*EPS)

  # v_obs has norm 1:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,v_obs),1.0,10*EPS)

  # v_obs is orthogonal to rho:
  assert_equal_eps(vector.inner_product(spacetime,chart,pars,x_obs,v_obs,rho),0.0,10*EPS)

  # rho has a positive r component:
  if r>1.0:
    assert_boolean(rho[1]>0.0,"observer's radial vector has a negative r component in Schwarzschild coordinates")

  # rho has vanishing angular components, observer is supposed to be radially infalling:
  assert_equal_eps(rho[2],0.0,10*EPS)
  assert_equal_eps(rho[3],0.0,10*EPS)
  assert_equal_eps(rho[4],0.0,10*EPS)

  # v_obs is oriented in the future-timelike direction:
  assert_boolean(transform.sch_is_in_future_light_cone(x_obs,v_obs),"v_obs is not future-oriented")
 
r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
test_obs_and_alpha(r)
r = 0.5 # test inside horizon
test_obs_and_alpha(r) # ... fails, possibly std observer is messed up for r<1?

test_riazuelo_deflection()
test_alpha_max_schwarzschild()
