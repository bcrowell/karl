#!/usr/bin/python3

# This file is basically just the test harness.

import sys
import numpy
numpy.seterr(all='raise')
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

import schwarzschild,util,io_util,sph_point,runge_kutta
from sph_point import SphPoint
from sph_vector import SphVector
from io_util import strcat,print_no_newline

def main():
  verbosity = 1
  do_test(verbosity,test_ks_sch_round_trip(verbosity))
  do_test(verbosity,test_ks_metric_against_sch_metric(verbosity))
  do_test(verbosity,test_ks_era(verbosity))
  do_test(verbosity,test_rotate_unit_sphere(verbosity))
  do_test(verbosity,test_create_sph_point(verbosity))
  do_test(verbosity,test_ks_christoffel_vs_raw_maxima(verbosity))
  do_test(verbosity,test_circular_orbit_dumb(verbosity))
  do_test(verbosity,test_geodesic_rk_simple(verbosity))
  simple=True
  do_test(verbosity,test_geodesic_rk_free_fall_from_rest(verbosity,simple))
  simple=False
  do_test(verbosity,test_geodesic_rk_free_fall_from_rest(verbosity,simple))
  do_test(verbosity,test_geodesic_rk_elliptical_period_fancy(verbosity))
  do_test(verbosity+1,test_ks_sch_transition_elliptical(verbosity+1))

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

# The Chistoffel symbols are analytic functions of the KS coordinates, so testing
# at one randomly chosen point has unit probability of detecting errors.
def test_ks_christoffel_vs_raw_maxima(verbosity):
  results = [True,""]
  theta = 0.789
  v = 0.123
  w = 0.456
  results = record_subtest(verbosity,results,subtest_ks_christoffel_vs_raw_maxima(verbosity,v,w,theta))
  return summarize_test(results,"test_ks_christoffel_vs_raw_maxima",verbosity)

def test_geodesic_rk_free_fall_from_rest(verbosity,simple):
  results = [True,""]
  # Choose initial and final r values such that if we're using the fancy version of the Runge-Kutta
  # routine, a transition will be triggered from Sch. to Kruskal.
  r1 = 4.0
  r2 = 2.0
  if not (r1>SphPoint.TRANSITION_MIN_R and r2<SphPoint.TRANSITION_MIN_R):
    return [True,"values of r1 and r2 are not above and below SphPoint.TRANSITION_MIN_R for testing"]
  results = record_subtest(verbosity,results,subtest_geodesic_rk_free_fall_from_rest(
                   verbosity,simple,100,r1,r2,SphPoint.SCHWARZSCHILD_CHART))
  if simple:
    results = record_subtest(verbosity,results,subtest_geodesic_rk_free_fall_from_rest(
                   verbosity,simple,100,r1,r2,SphPoint.KRUSKAL_VW_CHART))
  return summarize_test(results,"test_geodesic_rk_free_fall_from_rest",verbosity)

# Compute proper time for radial free fall from rest in Schwarzschild metric at
# Schwarzschild radius a to radius b and compare with exact equation (MTW p. 824, eq. 31.10).
# When we're in the weak-field region, Runge-Kutta is exact with n=1, because motion is quadratic in time.
# If setting the chart to Kruskal, keep in mind that large r will cause overflow errors.
def subtest_geodesic_rk_free_fall_from_rest(verbosity,simple,n,a,b,chart):
  info = ""
  ok = True
  t = 0.0
  theta = pi/2
  phi = 0.0
  # Start by constructing position and velocity in Schwarzschild chart, then switch both to Kruskal if needed.
  x = SphPoint(SphPoint.SCHWARZSCHILD,SphPoint.SCHWARZSCHILD_CHART,t,a,theta,phi)
  #x.debug_transitions = True
  tdot = 1.0/sqrt(1.0-1.0/a) # value of dt/dlambda such that the affine parameter lambda will be proper time
  v = SphVector(x,[tdot,0.0,0.0,0.0]) # derivative of coordinates with respect to proper time
  if chart==SphPoint.KRUSKAL_VW_CHART and simple:
    x.force_chart(chart) 
  if verbosity>=2: info += strcat(["initial point: chart=",x.chart,", x=",str(x),"\n"])
  eta = acos(2.0*b/a-1.0) # auxiliary parameter defined by MTW; eta=0 at start, and this is the value at b
  tau = 0.5*a**1.5*(eta+sin(eta))
  ndebug=0
  if verbosity>=2: ndebug = 10
  if simple:
    z = runge_kutta.geodesic_rk_simple(x,v,tau,tau/n,ndebug)
  else:
    z = runge_kutta.geodesic_rk       (x,v,tau,tau/n,ndebug,0,10,verbosity>=2) # ... ndebug_inner,ntrans,debug_transitions
  err = z[0]
  if err:
    print("error, "+z[1])
    exit(-1)
  final = z[2].absolute_schwarzschild()
  r = final[1]
  rel_err = (r-b)/b
  if verbosity>=2:
    info += strcat(["final tau=",tau,", r=",r,", expected r=",b," n=",n,", rel err=",rel_err,"\n"])
  eps = 1000.0/(n**4)
  if abs(rel_err)>eps:
    info += strcat(["relative discrepancy of ",rel_err," in final radius is greater than ",eps])
    ok = False
  return [ok,info]

# Test an elliptical orbit that causes transitions back and forth between Kruskal and Schwarzschild coordinates.
def test_ks_sch_transition_elliptical(verbosity):
  results = [True,""]
  r = 2.9
  if not (r<SphPoint.TRANSITION_MIN_R): return [False,"inappropriate perihelion radius"]
  a = 1.9 # elliptical
  f = 0.1 # fraction of the keplerian circular-orbit period
  method = 2 # makes it expressed as a fraction of the circular-orbit period
  simple = False
  n=1000
  results = record_subtest(verbosity,results,subtest_geodesic_rk_conserved(verbosity,n,r,a,f,method,simple))
  return summarize_test(results,"test_ks_sch_transition_elliptical",verbosity)

def test_geodesic_rk_simple(verbosity):
  results = [True,""]
  results = record_subtest(verbosity,results,subtest_geodesic_rk_simple_circular(verbosity))
  r = 1.0e8 # newtonian regime
  results = record_subtest(verbosity,results,subtest_geodesic_rk_elliptical_period(verbosity,100,r,1.0,True,0.0))
  results = record_subtest(verbosity,results,subtest_geodesic_rk_elliptical_period(verbosity,100,r,1.1,True,0.0))
  r = 10
  results = record_subtest(verbosity,results,subtest_geodesic_rk_conserved(verbosity,1000,r,1.1,0.3,1,True))
  return summarize_test(results,"test_geodesic_rk_simple",verbosity)

def test_geodesic_rk_elliptical_period_fancy(verbosity):
  results = [True,""]
  n = 100
  r = 1.0e8 # newtonian regime
  d = 0.0 # in equatorial plane

  a = 1.0 # circular orbit
  results = record_subtest(verbosity,results,subtest_geodesic_rk_elliptical_period(verbosity,n,r,a,False,d))

  a = 1.1 # elliptical orbit
  results = record_subtest(verbosity,results,subtest_geodesic_rk_elliptical_period(verbosity,n,r,a,False,d))

  d = 0.4 # inclined orbit, exercises other Christoffel symbols without rot90 transition
  results = record_subtest(verbosity,results,subtest_geodesic_rk_elliptical_period(verbosity,n,r,a,False,d))

  n = 1000 # with n=100, we get rather large errors in the polar orbit
  d = pi/2.0 # polar orbit, exercises rot90 transition
  results = record_subtest(verbosity,results,subtest_geodesic_rk_elliptical_period(verbosity,n,r,a,False,d))

  return summarize_test(results,"test_geodesic_rk_elliptical_period_fancy",verbosity)


# Only really tests two of the nonzero Christoffel symbols.
# The results are exact because C. symbols exactly cancel.
def subtest_geodesic_rk_simple_circular(verbosity):
  info = ""
  ok = True
  t = 0.0
  r = 10.0
  theta = pi/2
  phi = 0.0
  x = SphPoint(SphPoint.SCHWARZSCHILD,SphPoint.SCHWARZSCHILD_CHART,t,r,theta,phi)
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  v = SphVector(x,[1.0,0.0,0.0,v_phi]) # derivative of coordinates with respect to proper time
  if verbosity>=2: info += strcat(["initial point: chart=",x.chart,", x=",str(x),"\n"])
  period = 2.0*pi/v_phi
  n=10 # doesn't actually matter what n is, because it's exact
  ndebug=0
  if verbosity>=2: ndebug = 10
  z = runge_kutta.geodesic_rk_simple(x,v,period,period/n,ndebug)
  err = z[0]
  if err:
    print("error, "+z[1])
    exit(-1)
  final = z[2].absolute_schwarzschild()
  if abs(final[1]-r)>1.0e-12 or abs(final[2]-theta)>1.0e-12 or abs(sin(final[3])-sin(phi))>1.0e-12:
    info += strcat(["not back at starting position, final x=",str(x)])
    ok = False
  return [ok,info]

# Elliptical orbit.
# Start at perihelion. Make the initial velocity greater than the circular-orbit value by the factor a.
# Test against the Keplerian period. There is no point in testing with large n, because the errors
# become dominated by the Keplerian approximation.
# direction = angle about the x axis for the initial motion, defines plane of orbit
def subtest_geodesic_rk_elliptical_period(verbosity,n,r,a,simple,direction):
  info = ""
  ok = True
  t = 0.0
  theta = pi/2
  phi = 0.0
  x = SphPoint(SphPoint.SCHWARZSCHILD,SphPoint.SCHWARZSCHILD_CHART,t,r,theta,phi)
  if verbosity>=2: x.debug_transitions = True
  # Start by constructing parameters for a circular orbit, which we will later alter to be elliptical:
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  circular_period = 2.0*pi/v_phi
  # Increase velocity at perihelion to make orbit elliptical:
  v_phi = v_phi*a
  # Compute newtonian r_max:
  q=a**2/2.0-1.0
  r_max = r*((-1.0-sqrt(1.0+2.0*a**2*q))/(2*q))
  period = circular_period*((r+r_max)/(r+r))**1.5 # Kepler's law of periods
  v = SphVector(x,[1.0,0.0,v_phi*sin(direction),v_phi*cos(direction)]) 
       # ... derivative of coordinates with respect to proper time
  if verbosity>=2: info += strcat(["initial point: chart=",x.chart,", x=",str(x),"\n"])
  ndebug = 0
  if verbosity>=2: ndebug = 10
  if simple:
    z = runge_kutta.geodesic_rk_simple(x,v,period,period/n,ndebug)
  else:
    z = runge_kutta.geodesic_rk       (x,v,period,period/n,ndebug,0,10,verbosity>=2) # ... ndebug_inner,ntrans,debug_transitions
  err = z[0]
  if err:
    print("error, "+z[1])
    exit(-1)
  final = z[2].absolute_schwarzschild()
  eps = 100.0/r + 1000.0/(n**4) # first term is for error in Keplerian period, second for Runge-Kutta
  not_back = ( abs(final[1]/r-1.0)>eps or abs(final[2]-theta)>eps or abs(sin(final[3])-sin(phi))>eps )
  if not_back:
    info += strcat(["not back at starting position,\n  final (t,r,theta,phi)=",str(final),
                         "\n  errors in r,theta,phi=",
                         final[1]/r-1.0," ",final[2]-theta," ",sin(final[3])-sin(phi),"\n"])
    ok = False
  if verbosity>=2:
    info += strcat(["final point: chart=",x.chart,", x=",str(x),"\n"])
    info += strcat(["n=",n,", difference in sin phi=",abs(sin(final[3])-sin(phi)),"\n"])
    if ok:
      info += "passed"
    else:
      info += "failed"
  return [ok,info]

# Start at perihelion. Make the initial velocity greater than the circular-orbit value by the factor a.
# Test exactly conserved quantities. Go through a fraction f of the Keplerian estimate of the period,
# and check conserved quantities. Testing with large n makes sense, because these quantities are
# exactly relativistically conserved. However, for very large n we're dominated by rounding errors.
def subtest_geodesic_rk_conserved(verbosity,n,r,a,f,method,simple):
  info = ""
  ok = True
  t = 0.0
  theta = pi/2 # formula for angular momentum only works in the equatorial plane
  phi = 0.0 
  if method==1 and a>sqrt(2.0): raise RuntimeError('inappropriate value of a, not a bound orbit; use method=2 if you want this')
  x = SphPoint(SphPoint.SCHWARZSCHILD,SphPoint.SCHWARZSCHILD_CHART,t,r,theta,phi)
  x.debug_transitions = (verbosity>=2)
  # Start by constructing parameters for a circular orbit, which we will later alter to be elliptical:
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  circular_period = 2.0*pi/v_phi
  # Increase velocity at perihelion to make orbit elliptical:
  v_phi = v_phi*a
  if method==1:
    # Compute newtonian r_max:
    q=a**2/2.0-1.0
    r_max = r*((-1.0-sqrt(1.0+2.0*a**2*q))/(2*q))
    period = circular_period*((r+r_max)/(r+r))**1.5 # Kepler's law of periods
  else:
    # Method 2 is for use with unbound orbits or orbits that would be unbound in newtonian gravity.
    # Here we just use the keplerian period as a reasonable scale factor.
    period = circular_period
  max_lambda = f*period
  v = SphVector(x,[1.0,0.0,0.0,v_phi]) # derivative of coordinates with respect to proper time
  if not v.timelike():
    raise RuntimeError('velocity vector is not timelike')
  z = conserved_sch_stuff(x.get_raw_coords(),v.comp)
  l0 = z[0]; e0 = z[1]  
  if verbosity>=2: info += strcat(["initial point: chart=",x.chart,", x=",str(x),"\n"])
  if verbosity>=2: info += strcat(["initial values: L0=",("%8.5e" % l0),", E0=",("%8.5e" % e0),"\n"])
  ndebug = 0
  if verbosity>=2: ndebug = 10  
  if simple:
    z = runge_kutta.geodesic_rk_simple(x,v,max_lambda,max_lambda/n,ndebug)
  else:
    z = runge_kutta.geodesic_rk       (x,v,max_lambda,max_lambda/n,ndebug,0,10,verbosity>=2)
                                                     # ... ndebug_inner,ntrans,debug_transitions
  err = z[0]
  if err:
    print("error, "+z[1])
    exit(-1)
  final = z[2].absolute_schwarzschild()
  if verbosity>=2: info += strcat(["final point: chart=",x.chart,", x=",str(x),"\n"])
  z = conserved_sch_stuff(x.get_raw_coords(),v.comp)
  l = z[0]; e = z[1]  
  l_err = (l-l0)/l0
  e_err = (e-e0)/e0
  if verbosity>=2:
    info += strcat(["final values: L=",("%8.5e" % l),", E=",("%8.5e" % e),"\n"])
    info += strcat(["L err=",("%8.5e" % l_err ),", E err=",("%8.5e" % e_err )])
  eps = 1000.0/(n**4)
  if abs(l_err)>eps or abs(e_err)>eps:
    info += strcat(["conserved quantities not conserved, relative errors=",l_err," ",e_err])
    ok = False
  return [ok,info]

# Helper function for subtest_geodesic_rk_simple_conserved.
# https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Conserved_momenta
def conserved_sch_stuff(coords,v):
  r = coords[1]
  l = r*r*v[3]
  e = (1-1/r)*v[0]
  return [l,e]

def test_circular_orbit_dumb(verbosity):
  results = [True,""]
  results = record_subtest(verbosity,results,subtest_circular_orbit_dumb(verbosity))
  return summarize_test(results,"test_circular_orbit_dumb",verbosity)

# Dumb, low-tech test of whether we seem to get a circular orbit when we should.
# Doesn't use Runge-Kutta.
# Only really tests two of the nonzero Christoffel symbols.
# Despite the naive method for solving the ODEs, the results are exact because C. symbols exactly cancel.
def subtest_circular_orbit_dumb(verbosity):
  info = ""
  ok = True
  t = 0.0
  r = 10.0
  theta = pi/2
  phi = 0.0
  x = SphPoint(SphPoint.SCHWARZSCHILD,SphPoint.SCHWARZSCHILD_CHART,t,r,theta,phi)
  v_phi = 1/sqrt(2.0*r*r*r) # exact condition for circular orbit in Sch., if v_t=1.
  v = SphVector(x,[1.0,0.0,0.0,v_phi]) # derivative of coordinates with respect to proper time
  if verbosity>=2: info += strcat(["initial point: chart=",x.chart,", x=",str(x),"\n"])
  period = 2.0*pi/v_phi
  n = 10 # number of iterations; doesn't need to be big because solution is exact
  eps = 1.0/n
  dlambda = eps*period
  if verbosity>=2: info += strcat(["n=",n,", period=",period," dlambda=",dlambda," v=",v,"\n"])
  for iter in range(0,n):
    theta = x.theta
    ch = schwarzschild.sch_christoffel_sch(x.t,x.r,sin(theta),cos(theta))
    a = [0.0,0.0,0.0,0.0] # second derivative of x^i with respect to lambda
    for i in range(0, 4):
      for j in range(0, 4):
        for k in range(0, 4):
          a[i] = a[i] - ch[j][k][i]*v.comp[j]*v.comp[k]
    for i in range(0, 4):
      v.comp[i] += dlambda*a[i]
    x.add(dlambda,v)
  if verbosity>=2: info += strcat(["chart=",x.chart,", x=",str(x)])
  final = x.absolute_schwarzschild()
  if abs(final[1]-r)>1.0e-12 or abs(final[2]-theta)>1.0e-12 or abs(sin(final[3])-sin(phi))>1.0e-12:
    info += strcat(["not back at starting position, final x=",str(x)])
    ok = False
  return [ok,info]

def subtest_ks_christoffel_vs_raw_maxima(verbosity,v,w,theta):
  info = ""
  aux = schwarzschild.sch_aux_ks(v,w)
  r = aux[1]
  b = aux[2]
  ch = schwarzschild.sch_christoffel_ks(v,w,sin(theta),cos(theta),r,b)
  ch_raw = schwarzschild.sch_christoffel_ks_raw_maxima(v,w,theta,r)
  ok = True
  for i in range(0, 4):
    for j in range(0, 4):
      for k in range(0, 4):
        if abs(ch[i][j][k]-ch_raw[i][j][k])>1.0e-6:
          ok = False
          ratio = ch[i][j][k]/ch_raw[i][j][k]
          info += strcat(["disagreement for i=",i,", j=",j,", k=",k," raw=",ch_raw[i][j][k],", cooked=",ch[i][j][k],", ratio=",ratio,"\n"])
  return [ok,info]

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

if __name__ == '__main__':
  main()

###################################################################
