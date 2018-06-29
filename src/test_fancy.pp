#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector

verbosity=1

#--------------------------------------------------------------------------------------------------

def test_hitting_singularity():
  """
  Start from inside the horizon, so that Schwarzschild parameters work throughout. Simulate the
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
  v = vector.normalize(spacetime,chart,x,[dt,dr,0.0,0.0,0.0]) # normalize so that affine param is proper time
  tau_theory = (2.0/3.0)*r**1.5
  # ... theoretical time to hit the singularity for a trajectory with E=1
  #     https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Orbits_of_test_particles
  tau_max = 1.1*tau_theory
  n = 1000
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  opt = {'lambda_max':tau_max,'dlambda':tau_max/n,'ndebug':ndebug}
  err,final_x,final_v,final_a,final_lambda,info  = fancy.trajectory_schwarzschild(spacetime,chart,x,v,opt)
  assert_rel_equal_eps(final_lambda,tau_theory,1.0e-5)
  if err & RK_ERR:
    THROW('error: '+info['message'])


def main():
  test_hitting_singularity()

main()
