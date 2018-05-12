#!/usr/bin/python

# apt-get install python-scipy python-numpy

import sys
import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign
from scipy.special import lambertw

def main():
  verbosity = 1
  do_test(verbosity,test_ks_sch_round_trip(verbosity))
  do_test(verbosity,test_ks_metric_against_sch_metric(verbosity))
  do_test(verbosity,test_ks_era(verbosity))
  do_test(verbosity,test_rotate_unit_sphere(verbosity+1))

def do_test(verbosity,results):
  ok = results[0]
  info = results[1]
  if not ok :
    print("error in test!!!")
  print_no_newline(info)
  if not ok :
    print("-------------> exiting due to error")
    exit(-1)

def print_no_newline(s):
  sys.stdout.write(s)
  sys.stdout.flush()

def test_rotate_unit_sphere(verbosity):
  results = [True,""]
  theta = 0.784
  phi = 0.083 # ... needs to be small because:
  # The following round-trip conversions would not work if the input values of phi were not in the canonical range of [0,2pi).
  results = record_subtest(verbosity,results,subtest_rotate_unit_sphere_round_trip(verbosity,theta,phi))
  results = record_subtest(verbosity,results,subtest_rotate_unit_sphere_round_trip(verbosity,theta,phi+0.5*pi))
  results = record_subtest(verbosity,results,subtest_rotate_unit_sphere_round_trip(verbosity,theta,phi+pi))
  results = record_subtest(verbosity,results,subtest_rotate_unit_sphere_round_trip(verbosity,theta,phi+1.5*pi))
  return summarize_test(results,"test_rotate_unit_sphere",verbosity)
  

def test_ks_metric_against_sch_metric(verbosity):
  results = [True,""]
  theta = 0.784
  phi = 0.183
  results = record_subtest(verbosity,results,subtest_ks_metric_against_sch_metric(verbosity,0.1,-0.1,theta,phi))
                 # ... region I
  results = record_subtest(verbosity,results,subtest_ks_metric_against_sch_metric(verbosity,0.1,0.1,theta,phi))
                 # ... region II
  results = record_subtest(verbosity,results,subtest_ks_metric_against_sch_metric(verbosity,-0.1,0.1,theta,phi))
                 # ... region III
  results = record_subtest(verbosity,results,subtest_ks_metric_against_sch_metric(verbosity,-0.1,-0.1,theta,phi))
                 # ... region IV
  return summarize_test(results,"test_ks_metric_against_sch_metric",verbosity)

def test_ks_sch_round_trip(verbosity):
  results = [True,""]
  theta = 0.784
  phi = 0.183
  results = record_subtest(verbosity,results,subtest_ks_sch_round_trip(verbosity,0.1,-0.1,theta,phi)) # region I
  results = record_subtest(verbosity,results,subtest_ks_sch_round_trip(verbosity,0.1,0.1,theta,phi)) # region II
  results = record_subtest(verbosity,results,subtest_ks_sch_round_trip(verbosity,-0.1,0.1,theta,phi)) # region III
  results = record_subtest(verbosity,results,subtest_ks_sch_round_trip(verbosity,-0.1,-0.1,theta,phi)) # region IV
  return summarize_test(results,"test_ks_sch_round_trip",verbosity)

def test_ks_era(verbosity):
  results = [True,""]
  results = record_subtest(verbosity,results,subtest_ks_era(verbosity,0.1,-0.17)) # region I
  results = record_subtest(verbosity,results,subtest_ks_era(verbosity,0.1,0.17)) # region II
  results = record_subtest(verbosity,results,subtest_ks_era(verbosity,-0.1,0.17)) # region III
  results = record_subtest(verbosity,results,subtest_ks_era(verbosity,-0.1,-0.17)) # region IV
  return summarize_test(results,"test_ks_era",verbosity)

def subtest_ks_era(verbosity,v,w):
  info = ""
  t=3.7 # offset
  tx = ks_tx(v,w)
  old_sch = ks_to_sch(v,w)
  old_t = t+old_sch[0]
  ks = force_ks_era([t,v,w],tx,ks_to_region(v,w))
  new_sch = ks_to_sch(ks[1],ks[2])
  new_t = ks[0]+new_sch[0] # new_sch[0] should be 0, and we check that below
  if verbosity>=2: info += strcat(["old (t,V,W)=",t,",",v,",",w,"), new (t,V,W)=",ks[0],",",ks[1],",",ks[2],"), old_t=",old_t,",  new_t=",new_t])
  ok = (new_sch[0]==0.0) and (abs(new_t-old_t)<1.0e-8)
  return [ok,info]

def subtest_ks_metric_against_sch_metric(verbosity,v,w,theta,phi):
  info = ""
  sch = ks_to_sch(v,w)
  t = sch[0]
  r = sch[1]
  dv = 1.0e-5
  dw = 1.37e-5
  v2 = v+dv
  w2 = w+dw
  sch2 = ks_to_sch(v2,w2)
  t2 = sch2[0]
  r2 = sch2[1]
  dt = t2-t
  dr = r2-r
  g = sch_metric_ks(v,w,theta,phi)
  interval_ks = (g[0][1]+g[1][0])*dv*dw
  g = sch_metric_sch(r,theta)
  interval_sch = g[0][0]*dt**2+g[1][1]*dr**2
  if verbosity>=2: info += strcat(["V=",v,", W=",w," r=",r," t=",t,", int_ks=",interval_ks," int_sch=",interval_sch])
  ok = abs(interval_ks-interval_sch)<1.0e-13
  return [ok,info]

def subtest_ks_sch_round_trip(verbosity,v,w,theta,phi):
  info = ""
  sch = ks_to_sch(v,w)
  t = sch[0]
  r = sch[1]
  sigma = ks_to_sigma(v,w)
  # test round-trip transformation:
  rt_ks = sch_to_ks(t,r,sigma)
  rt_v = rt_ks[0]
  rt_w = rt_ks[1]
  if verbosity>=2: info += strcat(["V=",v,", W=",w,"; after round-trip transformation, V=",rt_v,", W=",rt_w])
  ok = (abs(v-rt_v)<1.0e-6) and (abs(w-rt_w)<1.0e-6)
  return [ok,info]

def subtest_rotate_unit_sphere_round_trip(verbosity,theta,phi):
  info = ""
  angles0 = [theta,phi]
  angles1 = rotate_unit_sphere(angles0,1.0)
  angles2 = rotate_unit_sphere(angles1,-1.0)
  theta2 = angles2[0]
  phi2 = angles2[1]
  if verbosity>=2: info += strcat([angles0," -> ",angles1," -> ",angles2])
  ok = (abs(theta-theta2)<1.0e-6) and (abs(phi-phi2)<1.0e-6)
  return [ok,info]


def summarize_test(results,name,verbosity):
  ok = results[0]
  if verbosity>=1 or not ok:
    if ok:
      results[1] += ("test "+name+" passed\n")
    else:
      results[1] += ("test "+name+" failed!!!!!!!!!!!!!!!\n")
    if verbosity>=2: results[1] += "============================\n"
  return results


def record_subtest(verbosity,results,subtest_results):
  ok = results[0]
  info = results[1]
  ok = (ok and subtest_results[0])
  info += (subtest_results[1]+"\n")
  if not ok:
    info += " ***FAILED***"
  #if verbosity>=2 or not ok: info += "\n"
  return [ok,info]  
    
# Define coordinate charts for different "eras," i.e., ranges of Schwarzschild t.
# Representation is (t,V,W), where t represents a time translation to be applied
# to the point (V,W), and we try to keep (V,W) such that it is fairly close to t=0.
# This prevents (V,W) from getting too big to represent in floating point. (They grow 
# exponentially with |t|.)
# This routine checks whether we've drifted out of the canonical chart, and if we have,
# fixes it.
def ks_era_in_range(ks):
  t = ks[0]
  v = ks[1]
  w = ks[2]
  tx = ks_tx(v,w)
  if abs(tx)<5.0:
    return ks # no change needed, already in canonical chart
  else:
    return force_ks_era(ks,tx,ks_to_region(v,w))

# Force (t,V,W) into a canonical form where (V,W) corresponds to a Schwarzschild t=0,
# i.e., the t information is all in the t parameter.
def force_ks_era(ks,tx,region):
  t = ks[0]
  v = ks[1]
  w = ks[2]
  t = t + 2.0*arctanh(tx)
  is_exterior = ks_is_exterior(region)
  if is_exterior:
    root_rho = sqrt(-v*w)
  else:
    root_rho = sqrt(v*w)
  v=root_rho*sign(v)
  w=root_rho*sign(w)
  return [t,v,w]

def ks_tx(v,w):
  region = ks_to_region(v,w)
  is_exterior = ks_is_exterior(region)
  # calculate T/X if exterior or X/T if interior:
  if is_exterior:
    return (v+w)/(v-w)
  else:
    return (v-w)/(v+w)


def ks_to_region(v,w):
  if v>=0 and w<=0: return 1 # exterior
  if v>=0 and w>=0: return 2 # interior
  if v<=0 and w>=0: return 3 # exterior, other copy of Minkowski space
  return 4 # white hole

def ks_is_exterior(region):
  return (region==1 or region==3)

def ks_to_sigma(v,w):
  return 1.0 if v>0 else -1.0  

# Take null Kruskal-Szekeres coordinates (V,W,theta,phi) and convert them 
# to Schwarzschild coordinates (t,r,theta,phi).
# The V,W coordinates are equivalent to Hawking and Ellis's (v'/sqrt2,w'/sqrt2).
# This will inevitably misbehave at the horizon, since KS coordinate are eliminating
# a coordinate singularity.
# This is my own way of expressing this transformation, may contain mistakes, needs to be tested.
def ks_to_sch(v,w):
  rho = -v*w
  l = lambert_w(rho/math.e)
  r = 1+l
  ks_x = (v-w)/2.0
  ks_t = (v+w)/2.0
  if rho>0.0:
    t = 2.0*arctanh(ks_t/ks_x) # exterior
  else:
    t = 2.0*arctanh(ks_x/ks_t) # interior
  return [t,r]

# This is guaranteed to misbehave at the horizon.
# sigma=+1 for regions I and II, -1 for III and IV
def sch_to_ks(t,r,sigma):
  t2 = 0.5*t
  if r>1.0:
    region=1 if sigma>1.0 else 3
    sc = sinh(t2)
    cs = cosh(t2)
  else:
    region=2 if sigma>1.0 else 4
    sc = cosh(t2)
    cs = sinh(t2)
  # My formulation based on MTW p. 827:
  a = sigma*sqrt(abs(r-1.0))*exp(r/2.0)
  ks_t = a*sc
  ks_x = a*cs
  v = ks_t+ks_x
  w = ks_t-ks_x
  return [v,w]

# For the Schwarzschild spacetime, compute the metric, given Kruskal-Szekeres null coordinates.
# Metric is in lower-index form, in +--- signature, with coordinates (V,W,theta,phi).
# The V,W coordinates are equivalent to Hawking and Ellis's (v'/sqrt2,w'/sqrt2).
def sch_metric_ks(v,w,theta,phi):
  g = [[0 for i in range(4)] for j in range(4)]
  rho = -v*w
  l = lambert_w(rho/math.e)
  r = 1+l
  b = 4.0*l/((1+l)*rho)
  g[0][1] = 0.5*b  
  g[1][0] = 0.5*b  
  r2 = r**2
  g[2][2] = -r2
  g[3][3] = -r2*sin(theta)**2
  return g

# For the Schwarzschild spacetime, compute the metric, given Schwarzschild coordinates.
# Metric is in lower-index form, in +--- signature, with coordinates (t,r,theta,phi).
# The mass is assumed to be 1/2, so that r is in units of the Schwarzschild radius.
# Angles in radians.
# Metric taken from https://en.wikipedia.org/wiki/Schwarzschild_metric .
def sch_metric_sch(r,theta):
  g = [[0 for i in range(4)] for j in range(4)]
  a = 1.0-1.0/r
  r2 = r**2
  g[0][0] = a
  g[1][1] = -1.0/a
  g[2][2] = -r2
  g[3][3] = -r2*sin(theta)**2
  return g

# If direction=+1, rotate a point on the unit sphere 90 degrees about
# the x axis in the direction that is right-handed with respect to
# positive x, so that, e.g., (theta,phi)=(0,0) -> (90,-90).
# If direction=-1, rotate the opposite direction, which is the inverse
# transformation.
def rotate_unit_sphere(angles,direction):
  theta = angles[0]
  phi = angles[1]
  s = sin(theta)
  theta2 = acos(-direction*s*sin(phi))
  phi2 = atan2(direction*cos(theta),s*cos(phi))
  if phi2<0.0: phi2 += 2.0*pi # python's atan2 returns values in (-pi,pi), but I want (0,2pi)
  return [theta2,phi2]

###################################################################

def lambert_w(x):
  return lambertw(x).real

###################################################################

# l is an array of objects which may not be strings
def strcat(l):
  return ''.join(map(str, l))

###################################################################

main()

###################################################################
