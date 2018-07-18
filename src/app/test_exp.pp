#!/usr/bin/python3

#####################
#  testing experimental stuff
#####################

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild
#if "LANG" eq "python"
import sys,os
#endif

verbosity=0

csv=TRUE
csv_file = 'a.csv'

#--------------------------------------------------------------------------------------------------

def test_err(r,chart,tau,e,l,n):
  """
  Start from inside the horizon. Simulate the
  motion of a massive test particle with energy e and angular momentum l.
  Return the [r,j] reached in proper time tau.
  """
  spacetime = SP_SCH
  # First calculate the initial conditions in Schwarzschild coordinates:
  x = [0.0,r,1.0,0.0,0.0]
  aa = 1-1/r
  tdot = e/aa
  phidot = l/(r*r)
  rdot = -sqrt(e*e-aa*(l*l/(r*r)+1))
  v = vector.normalize(spacetime,CH_SCH,{},x,[tdot,rdot,0.0,phidot,0.0])
  if chart!=CH_SCH:
    x2 = transform.transform_point(x,spacetime,CH_SCH,{},chart)
    v2 = transform.transform_vector(v,x,spacetime,CH_SCH,{},chart)
    x = x2
    v = v2
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  opt = {'lambda_max':tau,'dlambda':tau/n,'ndebug':ndebug}
  err,final_x,final_v,final_a,final_lambda,info  = \
              runge_kutta.trajectory_simple(SP_SCH,chart,{},x,v,opt)
  if err!=0:
    THROW("error: "+err)
  if verbosity>=2:
    PRINT("final_x=",final_x)
  if chart!=CH_SCH:
    x2 = transform.transform_point(final_x,spacetime,chart,{},CH_SCH)
    r=x2[1]
  else:
    r=final_x[1]  
  return [r,final_x[3]]


def debug_function(iter,lam,dlambda,x,v,name):
#if "LANG" eq "python"
  if csv:
    with open(csv_file, 'a') as f:
      f.write(io_util.fl_n_decimals(log(1.0e-30+abs(0.317303449000164-lam)),16)+","+io_util.fl_n_decimals(log(x[1]),16)+"\n")
#endif
  PRINT(name," lam=",io_util.fl(lam),", i=",iter,", dlam=",io_util.fl(dlambda),\
                      " x[0]=",io_util.fl_n_decimals(x[0],2), \
                      " x[1]=",io_util.fl_n_decimals(x[1],2),\
                      " v[0]=",io_util.fl_n_decimals(v[0],2), \
                      " v[1]=",io_util.fl_n_decimals(v[1],5))

verbosity=0

def main():

  EXIT(0) # don't run this as part of the routine test suite

  r = 1.0e-2
  e = 2.0
  l = 0.3
  n0 = 10
  k = 30.0 # steps per decade
  decades = 3.0
  niter = int(decades*k)
  if l!=0:
    tau_bound = (2.0/5.0)*(1.0/l)*r**(5.0/2.0)
  else:
    tau_bound = (2.0/3.0)*r**(3.0/2.0)
  tau=0.1*tau_bound
  est = 9.5891722323790222e-03
  n = 17
  for i in range(1000):
    rf,j = test_err(r,CH_SCH,tau,e,l,n)
    err = rf-est
    if False:
      PRINT(n,",",log(n)/log(2.0),",",io_util.fl(log(abs(err)+1.0e-30)/log(10.0)),",",io_util.fl_n_decimals(rf,16),",", \
                    io_util.fl_n_decimals(j,16))
      sys.stdout.flush()
  print("done")

main()
