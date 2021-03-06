#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import kruskal,transform,runge_kutta,angular,vector,lambert_w_stuff

#if "LANG" eq "python"
import ctypes
import c_libs
#endif

def test_round_trip_ksk(a,b):
  # We don't expect this to work if we're in regions III or IV, which can't be represented in S coords.
  t,r,mu = kruskal.aux(a,b)
  a2,b2 = transform.schwarzschild_to_kruskal(t,r)
  if verbosity>=3:
    PRINT("ksk: a=",a,", b=",b,", t=",t,", r=",r,", a'=",a2,", b'=",b2)
  test.assert_rel_equal_eps(a,a2,10*EPS)
  test.assert_rel_equal_eps(b,b2,10*EPS)

def test_round_trip_sks(t,r):
  a,b = transform.schwarzschild_to_kruskal(t,r)
  t2,r2,mu = kruskal.aux(a,b)
  if verbosity>=3:
    PRINT("sks: t=",t,", r=",r,", a=",a,", b=",b,", t2=",t2,", r2=",r2)
  test.assert_rel_equal_eps(t,t2,10*EPS)
  test.assert_rel_equal_eps(r,r2,10*EPS)

def simple_free_fall():
  # Free fall from rest in a semi-newtonian region, where all the math can be evaluated without
  # special massaging to avoid overflows.
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
  err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(SP_SCH,CH_AKS,{},x,v,opt)
  if err & RK_ERR:
    THROW('error: '+info['message'])
  tf,rf = transform.kruskal_to_schwarzschild(final_x[0],final_x[1])
  delta_r = r0-rf
  m=0.5 # coordinates are such that mass=1/2
  g=m/r0**2 # approximate as constant accel
  delta_r_newtonian = 0.5*g*lambda_max**2
  test.assert_rel_equal_eps(delta_r,delta_r_newtonian,0.02)

def test_christoffel_raw_vs_massaged(a,b):
  p = [a,b,0.8,0.6,0.0]
  ch = kruskal.christoffel_raw_maxima_output(p)
  ch2 = kruskal.christoffel_massaged_maxima_output(p)
  for i in range(5):
    for j in range(5):
      for k in range(5):
        #PRINT("i=",i,", j=",j,", k=",k,", raw=",ch[i][j][k],", massaged=",ch2[i][j][k])
        test.assert_rel_equal_eps(ch[i][j][k],ch2[i][j][k],10*EPS)

def test_motion_kruskal_vs_schwarzschild(t0,r0,flip,theta,phi,v,duration):
  # initial point = t0,r0,theta,phi
  # initial velocity = v, expressed in (t,r,i,j,k) component form
  # If flip is true, then the Kruskal version goes to (-a,-b).
  # For convenience, v will automatically be normalized and made tangent to the 2-sphere.
  # The variable duration is the maximum affine parameter, in units of the newtonian
  # free-fall time required in order to hit the singularity from this radius.
  i,j,k = angular.theta_phi_to_ijk(theta,phi)
  lambda_max = r0**1.5*duration
  a0,b0 = transform.schwarzschild_to_kruskal(t0,r0)
  if flip:
    a0 = -a0
    b0 = -b0
  # --- Find initial position in both coordinate systems: ------------------------------------
  x0s = CLONE_ARRAY_OF_FLOATS(([t0,r0,i,j,k])) # initial point in Schwarzschild coordinates
  x0k = CLONE_ARRAY_OF_FLOATS(([a0,b0,i,j,k])) # ... in Kruskal
  # --- Find initial velocity vector in Schwarzschild coordinates: ---------------------------
  v0s = angular.make_tangent(x0s,v) # initial velocity in Schwarzschild coordinates
  v0s = vector.normalize(SP_SCH,CH_SCH,{},x0s,v0s)
  # --- Find initial velocity vector in Kruskal coordinates: ---------------------------------
  jac = transform.jacobian_schwarzschild_to_kruskal(t0,r0)
  v0k = CLONE_ARRAY_OF_FLOATS(v0s)
  v0k[0] = jac[0][0]*v0s[0]+jac[0][1]*v0s[1]
  v0k[1] = jac[1][0]*v0s[0]+jac[1][1]*v0s[1]
  v0k = vector.normalize(SP_SCH,CH_AKS,{},x0k,v0k)
  if flip:
    v0k[0] = -v0k[0]
    v0k[1] = -v0k[1]
  # ---
  if verbosity>=3:
    PRINT("----------------------")
    PRINT("x0s=",io_util.vector_to_str_n_decimals(x0s,8))
    PRINT("x0k=",io_util.vector_to_str_n_decimals(x0k,8))
    PRINT("v0s=",io_util.vector_to_str_n_decimals(v0s,8))
    PRINT("v0k=",io_util.vector_to_str_n_decimals(v0k,8))
  # ---
  for i in range(2): # 0 for Schwarzschild, 1 for Kruskal
    n = 100 # fails with n=1000, why?
    ndebug = 0
    if FALSE and verbosity>=3:
      PRINT("----------------------")
      ndebug=n/100
    opt = {'lambda_max':lambda_max,'dlambda':lambda_max/n,'ndebug':ndebug}
    if i==0:
      x=x0s
      v=v0s
      chart = CH_SCH
    else:
      x=x0k
      v=v0k
      chart = CH_AKS
    err,final_x,final_v,final_a,final_lambda,info  = runge_kutta.trajectory_simple(SP_SCH,chart,{},x,v,opt)
    if err & RK_ERR:
      THROW('error: '+info['message'])
    if i==0:
      xfs = CLONE_ARRAY_OF_FLOATS(final_x)
    else:
      xfk = CLONE_ARRAY_OF_FLOATS(final_x)
  # Convert final result of Kruskal calculations to Schwarzschild coordinates:
  tf,rf = transform.kruskal_to_schwarzschild(xfk[0],xfk[1])
  xfk[0] = tf
  xfk[1] = rf
  if verbosity>=3:
    PRINT("xfs=",io_util.vector_to_str_n_decimals(xfs,8))
    PRINT("xfk=",io_util.vector_to_str_n_decimals(xfk,8))
  eps = 1000.0/n**4
  test.assert_rel_equal_eps_vector(xfs,xfk,eps)

#if "LANG" eq "python"

def test_aux_python_vs_c(a,b):
  t,r,mu = kruskal.aux(a,b)
  t2 = ctypes.c_double(0.0)
  r2 = ctypes.c_double(0.0)
  mu2 = ctypes.c_double(0.0)
  c_libs.karl_c_lib.aux(ctypes.byref(t2),ctypes.byref(r2),ctypes.byref(mu2),ctypes.c_double(a),ctypes.c_double(b))
  if verbosity>=3:
    print("t=",t,", r=",r,", mu=",mu)
    print("t2=",t2,", r2=",r2,", mu2=",mu2)
  # The results are exactly equal, which is a little surprising, since the implementations of W are different.
  # But I tested by intentionally altering one function, and the test did catch it.
  test.assert_rel_equal_eps(t,t2.value,EPS)
  test.assert_rel_equal_eps(r,r2.value,EPS)
  test.assert_rel_equal_eps(mu,mu2.value,EPS)

# The following only works if j==k, for the reasons explained in two comments below.
def test_christoffel_python_vs_c(a,b,j,k,i):
  if j!=k:
    THROW('j!=k')
  p = numpy.asarray([a,b,0.8,0.6,0.0,  0.0,0.0,0.0,0.0,0.0])
  p[5+j] = 1.0
  p[5+k] = 1.0
  ch = kruskal.christoffel(p)
  ch1=ch[j][k][i] # if j!=k, then there are two such terms, and the effect is double this
  a = numpy.asarray([0.0,0.0,0.0,0.0,0.0])
  p = p.astype(numpy.float64)
  a = a.astype(numpy.float64)
  c_double_p = ctypes.POINTER(ctypes.c_double)
  c_libs.karl_c_lib.sum_christoffel_sch_aks(p.ctypes.data_as(c_double_p),a.ctypes.data_as(c_double_p))
  ch2 = -a[i] # if j!=k, then we have terms Gamma^i_jj and Gamma^i_kk contributing
  test.assert_rel_equal_eps(ch1,ch2,5*EPS)

#endif

#--------------------------------------------------------------------------
# run tests
#--------------------------------------------------------------------------

#if "LANG" eq "python"

a=0.1
b=0.234
test_aux_python_vs_c(a,b)
test_christoffel_python_vs_c(a,b,0,0,0)
test_christoffel_python_vs_c(a,b,1,1,0)
test_christoffel_python_vs_c(a,b,0,0,1)
test_christoffel_python_vs_c(a,b,1,1,1)

a=0.1
b=-0.234
test_aux_python_vs_c(a,b)
test_christoffel_python_vs_c(a,b,0,0,0)
test_christoffel_python_vs_c(a,b,1,1,0)
test_christoffel_python_vs_c(a,b,0,0,1)
test_christoffel_python_vs_c(a,b,1,1,1)

#endif

def test_christoffel_symmetry(ch):
  l = len(ch)
  for i in range(l):
    for j in range(l):
      for k in range(l):
        test.assert_equal(ch[i][j][k],ch[j][i][k])

i,j,k = angular.theta_phi_to_ijk(0.333,0.1776)
test_christoffel_symmetry(kruskal.christoffel([1.0,-3.3,i,j,k])) # check symmetry at a random point

# region II:
test_round_trip_ksk(0.5,0.5)
# ... Note that if we test this with coordinates like (a,b)=(epsilon,epsilon), we get
#     poor relative precision, because in Schwarzschild coordinates, r=1+epsilon.

# region I:
test_round_trip_sks(0.0,2.0)
test_round_trip_ksk(1.0,-1.0)
test_round_trip_sks(100.0,100.0)

test_christoffel_raw_vs_massaged(0.5,-0.5) # region I
test_christoffel_raw_vs_massaged(0.5,0.5) # region II
test_christoffel_raw_vs_massaged(-0.5,0.5) # region III
test_christoffel_raw_vs_massaged(-0.5,-0.5) # region IV

simple_free_fall()

#---------------------------------------------------------------
# a bunch of applications of test_motion_kruskal_vs_schwarzschild:
#---------------------------------------------------------------

t0 = 0.0
r0 = 2.0
theta = MATH_PI/2.0
phi = 0.0
v = [1.0,0.0,0.0,0.0,0.0] # initially at rest; this gets normalized later
duration = 0.2
test_motion_kruskal_vs_schwarzschild(t0,r0,FALSE,theta,phi,v,duration) # region I
test_motion_kruskal_vs_schwarzschild(t0,r0,TRUE,theta,phi,v,duration) # region III
r0 = 0.5
v = [0.0,-1.0,0.0,0.0,0.0]
test_motion_kruskal_vs_schwarzschild(t0,r0,FALSE,theta,phi,v,duration) # region II
test_motion_kruskal_vs_schwarzschild(t0,r0,TRUE,theta,phi,v,duration) # region IV

# A test with random values of everything, to exercise all Christoffel symbols:
r0 = 2.0
theta = 1.111
phi = 2.345
v = [1.0,0.01776,0.01066,0.01492,0.02001]
duration = 0.2
test_motion_kruskal_vs_schwarzschild(t0,r0,FALSE,theta,phi,v,duration)

# Large r and t, test whether we get any overflows:
t0 = 0.0
r0 = 1.0e8
theta = MATH_PI/2.0
phi = 0.0
q = 1/r0 # scale angular motion down by this amount to keep the motion from being FTL
#v = [1.0,0.01776,0.01066*q,0.01492*q,0.02001*q] # random initial motion
v = [1.0,0,0,0,0] # random initial motion
duration = 1.0e-10 # Kruskal coordinates are not well adapted to covering the motion of
#                    a nonrelativistic object for long times.
test_motion_kruskal_vs_schwarzschild(t0,r0,FALSE,theta,phi,v,duration)

test.done(verbosity,"kruskal")
