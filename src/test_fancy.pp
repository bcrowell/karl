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
verbosity=1

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
  v = vector.normalize(spacetime,chart,x,[dt,dr,0.0,0.0,0.0])
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
  opt = {'lambda_max':tau_max,'dlambda':tau_max/n,'ndebug':ndebug,'debug_function':debug_function}
  err,final_x,final_v,final_a,final_lambda,info,sigma  = fancy.trajectory_schwarzschild(spacetime,chart,x,v,opt,1)
  #err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(spacetime,chart,x,v,opt)
  if verbosity>=2:
    PRINT("final_x=",final_x,", final_lambda=",final_lambda,", tau_theory=",tau_theory,", err=",(final_lambda-tau_theory))
  eps = 1.0e-9 # independent of n, because in Keplerian coords u' is constant
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

def main():
#if "LANG" eq "python"
  if csv and os.path.isfile(csv_file):
    os.remove(csv_file)
#endif
  test_hitting_singularity()

main()
