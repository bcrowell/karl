#!/usr/bin/python3

#include "math.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"

import runge_kutta

#verbosity=2

# Free fall from rest, from 10 Schwarzschild radii. 
# Affine param is approximately but not exactly equal to proper time.
x = [0.0,10.0,1.0,0.0,0.0]
v = [1.0,0.0,0.0,0.0,0.0]
opt = {'lambda_max':3.0,'dlambda':0.1,'ndebug':0}

# Smoke test.
err,final_x,final_v,final_lambda,info  = runge_kutta.geodesic_simple(SP_SCH,CH_SCH,x,v,opt)
if err & RK_ERR: raise RuntimeError('error: '+info['message'])
# Angular coordinates shouldn't have changed:
test.assert_equal(x[2],final_x[2])
test.assert_equal(x[3],final_x[3])
test.assert_equal(x[4],final_x[4])

# Newtonian limit. Free fall from ~1 a.u.
r0 = 1.0e8
lambda_max = r0**1.5*0.01 # short compared to the time required to hit the singularity
x = [0.0,r0,1.0,0.0,0.0]
v = [1.0,0.0,0.0,0.0,0.0]
opt = {'lambda_max':lambda_max,'dlambda':lambda_max/100,'ndebug':0}
err,final_x,final_v,final_lambda,info  = runge_kutta.geodesic_simple(SP_SCH,CH_SCH,x,v,opt)
if err & RK_ERR: raise RuntimeError('error: '+info['message'])
rf = final_x[1]
delta_r = r0-rf
m=1/2 # coordinates are such that mass=1/2
g=m/r0**2 # approximate as constant accel
delta_r_newtonian = 0.5*g*lambda_max**2
rel_err = (delta_r-delta_r_newtonian)/delta_r_newtonian
expect_err = 2*delta_r_newtonian/r0 # expected rel error due to constant-accel approximation
if verbosity>=2:
  print("delta_r=",delta_r,", delta_r_newtonian=",delta_r_newtonian," rel err=",rel_err,", expected=",expect_err)
if abs(rel_err)>2*expect_err: raise RuntimeError('error in final r greater than expected')
v_newtonian = -sqrt(2*m*(1/rf-1/r0))
vf = final_v[1]
rel_err_v = (vf-v_newtonian)/v_newtonian
if verbosity>=2:
  print("vf=",vf,", vf_newtonian=",v_newtonian," rel err=",rel_err_v,", expected=",expect_err)
if abs(rel_err_v)>2*expect_err: raise RuntimeError('error in final v greater than expected')

done(verbosity,"runge_kutta")
