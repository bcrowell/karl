#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian

verbosity=1

#--------------------------------------------------------------------------------------------------

def test_hitting_singularity():
  """
  Start from inside the horizon, so that Keplerian coordinates work throughout. Simulate the
  motion of a massive test particle on a trajectory that we know in closed form. Find the proper
  time to hit the singularity, and compare against theory.
  """
  spacetime = SP_SCH
  # First calculate the initial conditions in Schwarzschild coordinates:
  r = 0.9
  x = [0.0,r,1.0,0.0,0.0]
  aa = 1-1/r
  dt = 1.0
  dr = -sqrt(aa**2-aa**3) # the value that gives E=1=mc^2, i.e., falling from rest at infinity
  v = vector.normalize(spacetime,CH_SCH,x,[dt,dr,0.0,0.0,0.0]) # normalize so that affine param is proper time
  tau_theory = (2.0/3.0)*r**1.5
  # ... theoretical time to hit the singularity for a trajectory with E=1
  #     https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Orbits_of_test_particles
  tau_max = 1.1*tau_theory
  # Convert Schwarzschild coordinates to Keplerian:
  x[1] = x[1]**(3.0/2.0)
  v[1] = v[1]*(3.0/2.0)*sqrt(r)
  chart = CH_KEP
  n = 100
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  tol = 1.0e-5 # is exact up to rounding errors, regardless of n
  opt = {'lambda_max':tau_max,'dlambda':tau_max/n,'ndebug':ndebug,'tol':tol}
  err,final_x,final_v,final_a,final_lambda,info,sigma  = fancy.trajectory_schwarzschild(spacetime,chart,x,v,opt,1)
  if verbosity>=2:
    PRINT("final_x=",final_x)
  assert_rel_equal_eps(final_lambda,tau_theory,tol) 
  if err & RK_ERR:
    THROW('error: '+info['message'])


def main():
  test_hitting_singularity()

main()
