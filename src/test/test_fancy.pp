#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

#if "LANG" eq "python"
import os,os.path
#endif
import runge_kutta,fancy,angular,vector

csv=FALSE
csv_file = 'a.csv'
verbosity=0

#--------------------------------------------------------------------------------------------------

def test_hitting_singularity():
  """
  Start from inside the horizon, so that Schwarzschild coordinates work throughout. Simulate the
  motion of a massive test particle on a trajectory that we know in closed form. Find the proper
  time to hit the singularity, and compare against theory.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  r = 0.9
  x = [0.0,r,1.0,0.0,0.0]
  aa = 1-1/r
  dt = 1.0
  dr = -sqrt(aa**2-aa**3) # the value that gives E=1=mc^2, i.e., falling from rest at infinity
  v = vector.normalize(spacetime,chart,{},x,[dt,dr,0.0,0.0,0.0])
  # ...Normalize so that affine param is proper time.
  #    Since E=(1-1/r)t', I think this means we're using t->-t inside the horizon...?
  tau_theory = (2.0/3.0)*r**1.5
  # ... theoretical time to hit the singularity for a trajectory with E=1
  #     https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Orbits_of_test_particles
  tau_max = 1.1*tau_theory
  n = 1000
  ndebug=0
  if verbosity>=3:
    ndebug=n/20
    PRINT("tau_theory=",tau_theory)
    PRINT("x=",x)
    PRINT("v=",v)
  tol = 1.0e-6
  opt = {'lambda_max':tau_max,'dlambda':tau_max/n,'ndebug':ndebug,'debug_function':debug_function,'tol':tol,\
             'sigma':1,'future_oriented':TRUE}
  err,final_x,final_v,final_a,final_lambda,info,sigma  = fancy.trajectory_schwarzschild(spacetime,chart,{},x,v,opt)
  if verbosity>=2:
    PRINT("in test_hitting_singularity(), final_x=",final_x,", final_lambda=",final_lambda,", tau_theory=",tau_theory,", err=",(final_lambda-tau_theory))
  eps = tol
  if err!=0 and err!=RK_INCOMPLETE:
    THROW('error: '+info['message'])
  assert_rel_equal_eps(final_lambda,tau_theory,eps)

def debug_function(iter,lam,dlambda,x,v,name):
#if "LANG" eq "python"
  if name=='KEP' and csv:
    with open(csv_file, 'a') as f:
      f.write(io_util.fl_n_decimals(lam,16)+","+io_util.vector_to_str_n_decimals(x,16)+","+io_util.vector_to_str_n_decimals(v,16)+"\n")
#endif
  if not csv:
    PRINT(name," lam=",io_util.fl(lam),", i=",iter,", dlam=",io_util.fl(dlambda),\
                      " x[0]=",io_util.fl_n_decimals(x[0],2), \
                      " x[1]=",io_util.fl_n_decimals(x[1],2),\
                      " v[0]=",io_util.fl_n_decimals(v[0],2), \
                      " v[1]=",io_util.fl_n_decimals(v[1],5))

#--------------------------------------------------------------------------------------------------

def circular_orbit_period(tol):
  """
  Period of a circular orbit, Schwarzschild coordinates.

  There is also a version of this using non-adaptive RK, in test_runge_kutta.
  """
  r = 3.0
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  x = [0.0,r,1.0,0.0,0.0]
  v = [1.0,0.0,0.0,v_phi,0.0]
  period = 2.0*MATH_PI/v_phi
  opt = {'lambda_max':period,'ndebug':0,'norm_final':FALSE,'tol':tol,\
             'sigma':1,'future_oriented':TRUE}
  err,final_x,final_v,final_a,final_lambda,info,sigma  = fancy.trajectory_schwarzschild(SP_SCH,CH_SCH,{},x,v,opt)
  if verbosity>=2:
    PRINT("final x=",io_util.vector_to_str_n_decimals(final_x,16))
  if err & RK_ERR:
    THROW('error: '+info['message'])
  # Failing one of the following tests typically means that we need to readjust k in fancy.choose_step_size_exterior().
  test.assert_equal_eps(x[2],final_x[2],tol)
  test.assert_equal_eps(x[3],final_x[3],tol)
  test.assert_equal_eps(x[4],final_x[4],tol)

def test_chart_transitions():
  """
  Test a null geodesic that originates outside the horizon and ends at the singularity. This should
  cause two transitions between coordinate charts. This is an example that caused a crash, is now
  a smoke test; if it doesn't crash, the test succeeded.
  """
  x = [0.0, 2.0, 1.0, 0.0, 0.0]
  v = [1.0, -0.21794494717703375, 0.0, 0.31819805153394637, 0.0]
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  lambda_max = 10.0 # terminates near 6.8
  tol = 1.0e-6
  ndebug=0
  if verbosity>=3:
    ndebug=10
  opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,'tol':tol}
  err,final_x,final_v,final_a,final_lambda,info,sigma  = \
                fancy.trajectory_schwarzschild(spacetime,chart,pars,x,v,opt)
  if err!=RK_INCOMPLETE:
    PRINT('error in test_chart_transitions, return status was not RK_INCOMPLETE')
    PRINT(strcat(["final lambda=",final_lambda,", err=",err]))
    THROW('error')
  test.assert_equal_eps(final_lambda,6.74657323085,0.01) # should actually be much more accurate than this?

def main():
#if "LANG" eq "python"
  if csv and os.path.isfile(csv_file):
    os.remove(csv_file)
#endif
  circular_orbit_period(1.0e-3)
  test_hitting_singularity()
  test_chart_transitions()

main()
