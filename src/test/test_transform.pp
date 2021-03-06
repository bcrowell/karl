#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "init.h"
#include "test.h"

#include "precision.h"

import kruskal,transform

def test_round_trip_ksk(a,b):
  # We don't expect this to work if we're in regions III or IV, which can't be represented in S coords.
  t,r,mu = kruskal.aux(a,b)
  a2,b2 = transform.schwarzschild_to_kruskal(t,r)
  if verbosity>=3:
    PRINT("ksk: a=",a,", b=",b,", t=",t,", r=",r,", a'=",a2,", b'=",b2)
  test.assert_rel_equal_eps(a,a2,5*EPS)
  test.assert_rel_equal_eps(b,b2,5*EPS)

def test_round_trip_sks(t,r):
  a,b = transform.schwarzschild_to_kruskal(t,r)
  t2,r2,mu = kruskal.aux(a,b)
  if verbosity>=3:
    PRINT("sks: t=",t,", r=",r,", a=",a,", b=",b,", t2=",t2,", r2=",r2)
  test.assert_rel_equal_eps(t,t2,5*EPS)  
  test.assert_rel_equal_eps(r,r2,5*EPS)  

def test_jacobian_schwarzschild_to_kruskal(t,r):
  # Test analytic Jacobian against numerical differentiation.
  jac = transform.jacobian_schwarzschild_to_kruskal(t,r)
  d = sqrt(EPS) # differences of 10^-8, so that we can get numerical approximations to derivatives to 8 decimals
  a0,b0 = transform.schwarzschild_to_kruskal(t,r)
  for i in range(2): # a,b
    for j in range(2): # t,r
      t1 = t
      r1 = r
      if j==0:
        t1 += d
      if j==1:
        r1 += d
      a1,b1 = transform.schwarzschild_to_kruskal(t1,r1)
      if i==0:
        dk=a1-a0
      if i==1:
        dk=b1-b0
      der = dk/d # numerical approximation to the partial derivative
      if verbosity>=3:
        PRINT("i=",i,", j=",j,", der=",der,", J=",jac[i][j])
      test.assert_rel_equal_eps(der,jac[i][j],10*d)

def test_jacobian_kruskal_to_schwarzschild(t,r):
  # Test analytic Jacobian against numerical differentiation.
  jac = transform.jacobian_kruskal_to_schwarzschild(t,r)
  d = sqrt(EPS) # differences of 10^-8, so that we can get numerical approximations to derivatives to 8 decimals
  a,b = transform.schwarzschild_to_kruskal(t,r)
  for i in range(2): # t,r
    for j in range(2): # a,b
      a1 = a
      b1 = b
      if j==0:
        a1 += d
      if j==1:
        b1 += d
      t1,r1 = transform.kruskal_to_schwarzschild(a1,b1)
      if i==0:
        ds=t1-t
      if i==1:
        ds=r1-r
      der = ds/d # numerical approximation to the partial derivative
      if verbosity>=3:
        PRINT("i=",i,", j=",j,", der=",der,", J=",jac[i][j])
      test.assert_rel_equal_eps(der,jac[i][j],10*d)

# region II:
test_round_trip_ksk(0.5,0.5)
# ... Note that if we test this with coordinates like (a,b)=(epsilon,epsilon), we get
#     poor relative precision, because in Schwarzschild coordinates, r=1+epsilon.

# region I:
test_round_trip_sks(0.0,2.0)
test_round_trip_ksk(1.0,-1.0)
test_round_trip_sks(100.0,100.0)

# In the following, the arguments are (t,r).
test_jacobian_schwarzschild_to_kruskal(0.111,2.0) # region I
test_jacobian_schwarzschild_to_kruskal(0.222,0.5) # region II
test_jacobian_kruskal_to_schwarzschild(0.111,2.0) # region I
test_jacobian_kruskal_to_schwarzschild(0.222,0.5) # region II


test.done(verbosity,"transform")
