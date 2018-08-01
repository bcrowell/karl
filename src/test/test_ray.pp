#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"
#include "spacetimes.h"

import test,vector,ray,star_properties

def test_obs_and_alpha(r):
  test_schwarzschild_standard_observer(r) 
  test_le_to_alpha_schwarzschild(r)

def test_riazuelo_deflection():
  le_over_max = 0.728
  alpha,beta = test_riazuelo_deflection_one(le_over_max)
  print("test_riazuelo_deflection(): alpha=",alpha,", beta=",beta)

def test_riazuelo_deflection_one(le_over_max):
  r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
  # Interpolating from his graph, he has lines crossing at alpha=beta=2.48+-0.1.
  # Working backwards, this requires about L/E=11.3, which is 0.728 of max:
  #     calc -e "r=15; aa=1-1/r; lemax=r/sqrtaa; le=11.3/lemax"
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  x_obs,v_obs,rho = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1-1/r
  le_max = r/sqrt(abs(aa)) # abs is so we still get a test at r<1, see above
  le = le_over_max*le_max
  in_n_out = 1
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
  print("r=",r,", le_over_max=",0.1,", alpha=",alpha)

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

  if chart==CH_SCH:
    # rho has a positive r component:
    if r>1.0:
      assert_boolean(rho[1]>0.0,"observer's radial vector has a negative r component in Schwarzschild coordinates")

    # rho has vanishing angular components, observer is supposed to be radially infalling:
    assert_equal_eps(rho[2],0.0,10*EPS)
    assert_equal_eps(rho[3],0.0,10*EPS)
    assert_equal_eps(rho[4],0.0,10*EPS)

    # v_obs is oriented in the future-timelike direction:
    if r>1.0:
      assert_boolean(v_obs[0]>0.0,"observer has r>1 and v_t<=0, i.e., velocity is not future-oriented")
    if r<1.0:
      assert_boolean(v_obs[1]<0.0,"observer has r<1 and v_t>=0, i.e., velocity is not future-oriented")

 
r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
test_obs_and_alpha(r)
r = 0.5 # test inside horizon
#test_obs_and_alpha(r) # ... fails, possibly std observer is messed up for r<1?

test_riazuelo_deflection()
