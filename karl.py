#!/usr/bin/python3

# This file is basically just the test harness.

# apt-get install python-scipy3 python-numpy3

import sys
import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

import schwarzschild,util,sph_point
from util import lambert_w
from sph_point import SphPoint

def main():
  verbosity = 1
  do_test(verbosity,test_ks_sch_round_trip(verbosity))
  do_test(verbosity,test_ks_metric_against_sch_metric(verbosity))
  do_test(verbosity,test_ks_era(verbosity))
  do_test(verbosity,test_rotate_unit_sphere(verbosity))
  do_test(verbosity+1,test_create_sph_point(verbosity+1))

def do_test(verbosity,results):
  ok = results[0]
  info = results[1]
  if not ok :
    print("error in test!!!")
  print_no_newline(info)
  if not ok :
    print("-------------> exiting due to error")
    exit(-1)

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

def test_create_sph_point(verbosity):  
  results = [True,""]
  theta = 0.7
  phi = 0.1
  r = 10.0 # a value that will not force conversion to Kruskal chart
  t = 0.0
  results = record_subtest(verbosity,results,subtest_create_sph_point(verbosity,t,r,theta,phi))
  theta = 0.0 # test at coordinate singularity, so a rotation happens
  results = record_subtest(verbosity,results,subtest_create_sph_point(verbosity,t,r,theta,phi))
  theta = pi
  results = record_subtest(verbosity,results,subtest_create_sph_point(verbosity,t,r,theta,phi))
  theta = 0.7
  r = 1.1 # a value that will force conversion to Kruskal chart
  results = record_subtest(verbosity,results,subtest_create_sph_point(verbosity,t,r,theta,phi))
  r = 10.0
  t = 10.0 # a value that will force the use of the time translation (ks_era_in_range)
  results = record_subtest(verbosity,results,subtest_create_sph_point(verbosity,t,r,theta,phi))
  return summarize_test(results,"test_create_sph_point",verbosity)

# All this does is exercise the code used in creating the object. Doesn't check if the results make sense.
# This is basically just a smoke test to make sure the code runs.
def subtest_create_sph_point(verbosity,t,r,theta,phi):
  info = ""
  p = SphPoint(SphPoint.SCHWARZSCHILD,SphPoint.SCHWARZSCHILD_CHART,t,r,theta,phi)
  k = p.absolute_kruskal()
  if verbosity>=2: info += strcat(["k=",k])
  ok = True
  return [ok,info]

def subtest_ks_era(verbosity,v,w):
  info = ""
  t=3.7 # offset
  tx = schwarzschild.ks_tx(v,w)
  old_sch = schwarzschild.ks_to_sch(v,w)
  old_t = t+old_sch[0]
  ks = schwarzschild.force_ks_era([t,v,w],tx,schwarzschild.ks_to_region(v,w))
  new_sch = schwarzschild.ks_to_sch(ks[1],ks[2])
  new_t = ks[0]+new_sch[0] # new_sch[0] should be 0, and we check that below
  if verbosity>=2: info += strcat(["old (t,V,W)=",t,",",v,",",w,"), new (t,V,W)=",ks[0],",",ks[1],",",ks[2],"), old_t=",old_t,",  new_t=",new_t])
  ok = (new_sch[0]==0.0) and (abs(new_t-old_t)<1.0e-8)
  return [ok,info]

def subtest_ks_metric_against_sch_metric(verbosity,v,w,theta,phi):
  info = ""
  sch = schwarzschild.ks_to_sch(v,w)
  t = sch[0]
  r = sch[1]
  dv = 1.0e-5
  dw = 1.37e-5
  v2 = v+dv
  w2 = w+dw
  sch2 = schwarzschild.ks_to_sch(v2,w2)
  t2 = sch2[0]
  r2 = sch2[1]
  dt = t2-t
  dr = r2-r
  aux = schwarzschild.sch_aux_ks(v,w)
  rho=aux[0]; r=aux[1]; b=aux[2]
  sin_theta = sin(theta)
  g = schwarzschild.sch_metric_ks(v,w,sin_theta,r,b)
  interval_ks = (g[0][1]+g[1][0])*dv*dw
  g = schwarzschild.sch_metric_sch(r,sin_theta)
  interval_sch = g[0][0]*dt**2+g[1][1]*dr**2
  if verbosity>=2: info += strcat(["V=",v,", W=",w," r=",r," t=",t,", int_ks=",interval_ks," int_sch=",interval_sch])
  ok = abs(interval_ks-interval_sch)<1.0e-13
  return [ok,info]

def subtest_ks_sch_round_trip(verbosity,v,w,theta,phi):
  info = ""
  sch = schwarzschild.ks_to_sch(v,w)
  t = sch[0]
  r = sch[1]
  sigma = schwarzschild.ks_to_sigma(v,w)
  # test round-trip transformation:
  rt_ks = schwarzschild.sch_to_ks(t,r,sigma)
  rt_v = rt_ks[0]
  rt_w = rt_ks[1]
  if verbosity>=2: info += strcat(["V=",v,", W=",w,"; after round-trip transformation, V=",rt_v,", W=",rt_w])
  ok = (abs(v-rt_v)<1.0e-6) and (abs(w-rt_w)<1.0e-6)
  return [ok,info]

def subtest_rotate_unit_sphere_round_trip(verbosity,theta,phi):
  info = ""
  angles0 = [theta,phi]
  angles1 = util.rotate_unit_sphere(angles0,1.0)
  angles2 = util.rotate_unit_sphere(angles1,-1.0)
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
  if subtest_results[1]!="": info += (subtest_results[1]+"\n") 
  if not ok:
    info += " ***FAILED***"
  return [ok,info]  

###################################################################

# l is an array of objects which may not be strings
def strcat(l):
  return ''.join(map(str, l))

def print_no_newline(s):
  sys.stdout.write(s)
  sys.stdout.flush()

###################################################################

if __name__ == '__main__':
  main()

###################################################################
