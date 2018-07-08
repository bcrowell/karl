#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"

import runge_kutta,angular,vector,transform

#--------------------------------------------------------------------------------------------------

def smoke_test():
  # Free fall from rest, from 10 Schwarzschild radii. 
  # Affine param is approximately but not exactly equal to proper time.
  x = [0.0,10.0,1.0,0.0,0.0]
  v = [1.0,0.0,0.0,0.0,0.0]
  ndebug=0
  if verbosity>=3:
    ndebug=1
  opt = {'lambda_max':100.0,'dlambda':10.0,'ndebug':ndebug}
  err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(SP_SCH,CH_SCH,x,v,opt)
  if err & RK_ERR:
    THROW('error: '+info['message'])
  # Angular coordinates shouldn't have changed:
  test.assert_equal(x[2],final_x[2])
  test.assert_equal(x[3],final_x[3])
  test.assert_equal(x[4],final_x[4])

#--------------------------------------------------------------------------------------------------

def test_force():
  # Hover at r Schwarzschild radii.
  # Affine param is approximately but not exactly equal to proper time.
  r = 1000.0
  x = [0.0,r,1.0,0.0,0.0]
  v = [1.0,0.0,0.0,0.0,0.0]
  ndebug=0
  if verbosity>=3:
    ndebug=1
  fmag = 0.5/(r*r) # Newtonian acceleration, m=1/2
  force = [0.0,fmag,0.0,0.0,0.0] # constant radial proper acceleration
  f = lambda lam,x,v: force #js f = function(lam,x,v) {return force;}
  opt = {'lambda_max':100.0,'dlambda':10.0,'ndebug':ndebug,\
         'force_acts':TRUE,'force_function':f,'force_chart':CH_SCH}
  err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(SP_SCH,CH_SCH,x,v,opt)
  if err & RK_ERR:
    THROW('error: '+info['message'])
  # We're hovering, so r shouldn't have changed:
  test.assert_rel_equal_eps(x[1],final_x[1],1.0e-8)
  # Angular coordinates shouldn't have changed:
  test.assert_equal(x[2],final_x[2])
  test.assert_equal(x[3],final_x[3])
  test.assert_equal(x[4],final_x[4])

#--------------------------------------------------------------------------------------------------

def simple_newtonian_free_fall():
  # Newtonian limit. Free fall from ~1 a.u.
  r0 = 1.0e8
  lambda_max = r0**1.5*0.01 # short compared to the time required to hit the singularity
  x = [0.0,r0,1.0,0.0,0.0]
  v = [1.0,0.0,0.0,0.0,0.0]
  n = 100
  opt = {'lambda_max':lambda_max,'dlambda':lambda_max/n,'ndebug':0}
  err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(SP_SCH,CH_SCH,x,v,opt)
  if err & RK_ERR:
    THROW('error: '+info['message'])
  rf = final_x[1]
  delta_r = r0-rf
  m=0.5 # coordinates are such that mass=1/2
  g=m/r0**2 # approximate as constant accel
  delta_r_newtonian = 0.5*g*lambda_max**2
  rel_err = (delta_r-delta_r_newtonian)/delta_r_newtonian
  expect_err = 2*delta_r_newtonian/r0 # expected rel error due to constant-accel approximation
  if verbosity>=2:
    PRINT("delta_r=",delta_r,", delta_r_newtonian=",delta_r_newtonian," rel err=",rel_err,", expected=",expect_err)
  test.assert_rel_equal_eps(delta_r,delta_r_newtonian,expect_err)
  if abs(rel_err)>2*expect_err:
    THROW('error in final r greater than expected')
  v_newtonian = -sqrt(2*m*(1/rf-1/r0))
  vf = final_v[1]
  rel_err_v = (vf-v_newtonian)/v_newtonian
  if verbosity>=2:
    PRINT("vf=",vf,", vf_newtonian=",v_newtonian," rel err=",rel_err_v,", expected=",expect_err)
  test.assert_rel_equal_eps(vf,v_newtonian,expect_err)
  if abs(rel_err_v)>2*expect_err:
    THROW('error in final v greater than expected')

#--------------------------------------------------------------------------------------------------

def circular_orbit_period():
  """
  Period of a circular orbit, Schwarzschild coordinates.
  """
  r = 3.0
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  x = [0.0,r,1.0,0.0,0.0]
  v = [1.0,0.0,0.0,v_phi,0.0]
  n = 100
  period = 2.0*MATH_PI/v_phi
  opt = {'lambda_max':period,'dlambda':period/n,'ndebug':0,'norm_final':FALSE}
  err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(SP_SCH,CH_SCH,x,v,opt)
  if verbosity>=2:
    PRINT("final x=",io_util.vector_to_str_n_decimals(final_x,16))
  if err & RK_ERR:
    THROW('error: '+info['message'])
  eps = 1.0e4/n**4
  test.assert_equal_eps(x[2],final_x[2],eps)
  test.assert_equal_eps(x[3],final_x[3],eps)
  test.assert_equal_eps(x[4],final_x[4],eps)

#--------------------------------------------------------------------------------------------------

def elliptical_orbit_period(r,a,direction,n,half_period):
  """
  Period of an elliptical orbit, Schwarzschild coordinates.

  Start at perihelion, r. Make the initial velocity greater than the circular-orbit value by the factor a.
  Test against the Keplerian period. There is no point in testing with large n, because the errors
  become dominated by the Keplerian approximation.
  direction = angle about the x axis for the initial motion, defines plane of orbit
  Runge-Kutta with n steps.
  This is intended to be used with very large r, so that the Keplerian approximation is good.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  x = [0.0,r,1.0,0.0,0.0]
  circular_period = 2.0*MATH_PI/v_phi
  # Increase velocity at perihelion to make orbit elliptical:
  v_phi = v_phi*a
  v = [1.0,0.0,0.0,v_phi*cos(direction),v_phi*sin(direction)]
  v = angular.make_tangent(x,v)
  v = vector.normalize(spacetime,chart,x,v)
  # Compute newtonian r_max:
  q=a**2/2.0-1.0
  r_max = r*((-1.0-sqrt(1.0+2.0*a**2*q))/(2*q))
  foo = 1.0+2.0*a**2*q
  period = circular_period*((r+r_max)/(r+r))**1.5 # Kepler's law of periods
  triggers = []
  lambda_max = period
  if half_period:
    lambda_max = 0.5*period*(1+2.0/n) # has to be longer than 0.5*period, or it doesn't test trigger
    triggers = [[-1.0, 6, 0, 0.5]] # trigger when rdot crosses zero from above
  #--
  ndebug=0
  if verbosity>=3:
    if half_period:
      PRINT("testing half-period")
    ndebug=n/10
  opt = {'lambda_max':lambda_max,'dlambda':lambda_max/n,'ndebug':ndebug,'norm_final':FALSE,'triggers':triggers}
  err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(spacetime,chart,x,v,opt)
  if err & RK_ERR:
    THROW('error: '+info['message'])
  if verbosity>=2:
    PRINT("final x=",io_util.vector_to_str_n_decimals(final_x,16))
  if half_period:
    eps = 10.0/n # won't be as accurate as RK, because we extrapolate linearly to find where v crosses zero
    lamx = -final_v[1]/final_a[1] # additional increment to lambda based on extrapolation to the trigger
    final_lambda = final_lambda+lamx
    final_t = final_x[0]+final_v[0]*lamx
    final_r = final_x[1] # any correction would be second order
    final_j = final_x[3]+final_v[3]*lamx
    if verbosity>=3:
      PRINT("final lam=",final_lambda,", t=",final_t,", r=",final_r,", j=",final_j)
    test.assert_rel_equal_eps(final_r,r_max,eps)
    test.assert_rel_equal_eps(final_t,0.5*period,eps)
    test.assert_equal_eps(final_j,0.0,eps)
  else:
    eps = 100.0/r + 10000.0/(n**4) # first term is for error in Keplerian period, second for Runge-Kutta
    test.assert_equal_eps(x[2],final_x[2],eps)
    test.assert_equal_eps(x[3],final_x[3],eps)
    test.assert_equal_eps(x[4],final_x[4],eps)


#--------------------------------------------------------------------------------------------------

def test_radial_null_geodesic(r,rdot,lam,chart,n):
  """
  Start a radial null geodesic at r with dr/dlambda=rdot.
  Return the r reached at affine parameter lam.
  Test against the closed-form solution.
  """
  spacetime = SP_SCH
  # First calculate the initial conditions in Schwarzschild coordinates:
  x = [0.0,r,1.0,0.0,0.0]
  aa = 1-1/r
  tdot = abs(rdot/aa) # Either sign is actually possible.
  v = [tdot,rdot,0.0,0.0,0.0]
  if chart!=CH_SCH:
    x2 = transform.transform_point(x,spacetime,CH_SCH,chart)
    v2 = transform.transform_vector(v,x,spacetime,CH_SCH,chart)
    x = x2
    v = v2
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  opt = {'lambda_max':lam,'dlambda':lam/n,'ndebug':ndebug}
  err,final_x,final_v,final_a,final_lambda,info  = \
              runge_kutta.trajectory_simple(SP_SCH,chart,x,v,opt)
  if err!=0:
    THROW("error: "+err)
  final_x = transform.transform_point(final_x,spacetime,chart,CH_SCH) # convert to Schwarzschild coords
  if verbosity>=2:
    PRINT("final_x (Schwarzschild)=",final_x)
  # Test against the closed-form solution r=r0+rdot*lam, t=const+-(r+ln|r-1|).
  rf = final_x[1]
  tf = final_x[0]
  eps = 1.0/n**4
  test.assert_equal_eps(rf,r+rdot*lam,eps)
  test.assert_equal_eps(abs(tf),abs((r+log(abs(r-1.0)))-(rf+log(abs(rf-1.0)))),eps)

#--------------------------------------------------------------------------------------------------

verbosity=1

def main():
  test_radial_null_geodesic(2.0,1.0,1.0,CH_SCH,100)
  smoke_test()
  simple_newtonian_free_fall()
  circular_orbit_period()
  #--
  r = 1.0e8
  a = 1.1
  direction = 0.0
  n = 100
  elliptical_orbit_period(r,a,direction,n,FALSE) # test period
  elliptical_orbit_period(r,a,direction,n,TRUE) # test half-period
  #--
  test_force()
  #--
  test.done(verbosity,"runge_kutta")

main()
