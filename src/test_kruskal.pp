#!/usr/bin/python3

#include "math.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import kruskal,transform,runge_kutta

def test_round_trip_ksk(a,b):
  # We don't expect this to work if we're in regions III or IV, which can't be represented in S coords.
  t,r,mu = kruskal.aux(a,b)
  a2,b2 = transform.schwarzschild_to_kruskal(t,r)
  if verbosity>=3: print("ksk: a=",a,", b=",b,", t=",t,", r=",r,", a'=",a2,", b'=",b2)
  test.assert_rel_equal(a,a2)
  test.assert_rel_equal(b,b2)

def test_round_trip_sks(t,r):
  a,b = transform.schwarzschild_to_kruskal(t,r)
  t2,r2,mu = kruskal.aux(a,b)
  if verbosity>=3: print("sks: t=",t,", r=",r,", a=",a,", b=",b,", t2=",t2,", r2=",r2)
  test.assert_rel_equal(t,t2)  
  test.assert_rel_equal(r,r2)  

def simple_free_fall():
  t0 = 0.0
  r0 = 300.0 # as big as possible without causing overflows
  lambda_max = r0**1.5*0.01 # short compared to the Newtonian time required to hit the singularity
  a,b = transform.schwarzschild_to_kruskal(t0,r0)
  x = [a,b,1.0,0.0,0.0]
  jac = transform.jacobian_schwarzschild_to_kruskal(t0,r0)
  tdot = 1.0 # approximate newtonian value of dt/dtau
  v = [jac[0][0]*tdot,jac[1][0]*tdot,0.0,0.0,0.0]
  n = 100
  opt = {'lambda_max':lambda_max,'dlambda':lambda_max/n,'ndebug':0}
  err,final_x,final_v,final_lambda,info  = runge_kutta.geodesic_simple(SP_SCH,CH_AKS,x,v,opt)
  if err & RK_ERR: raise RuntimeError('error: '+info['message'])
  tf,rf = transform.kruskal_to_schwarzschild(final_x[0],final_x[1])
  delta_r = r0-rf
  m=1/2 # coordinates are such that mass=1/2
  g=m/r0**2 # approximate as constant accel
  delta_r_newtonian = 0.5*g*lambda_max**2
  test.assert_rel_equal_eps(delta_r,delta_r_newtonian,0.02)


# region II:
test_round_trip_ksk(0.5,0.5)
# ... Note that if we test this with coordinates like (a,b)=(epsilon,epsilon), we get
#     poor relative precision, because in Schwarzschild coordinates, r=1+epsilon.

# region I:
test_round_trip_sks(0.0,2.0)
test_round_trip_ksk(1.0,-1.0)
test_round_trip_sks(100.0,100.0)

simple_free_fall()

done(verbosity,"kruskal")
