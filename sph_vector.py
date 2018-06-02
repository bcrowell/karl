import copy

import schwarzschild,util

import numpy
numpy.seterr(all='raise')
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

import sph_point
from sph_point import SphPoint

# Documentation for the math is in the file doc.tex, which can be
# compiled to pdf format by doing a "make doc." (Comments in the code do
# not document the math or the definitions of the variables.)

class SphVector:
  # This class represents a vector in the tangent space at a particular point.
  # The point is a member of the class SphPoint.
  # Although mathematically we think of a vector as being tied to the tangent space of one exact
  # point on the manifold, here we actually take the liberty of letting the point change
  # without doing anything to the vector, provided that we didn't change to a different chart.
  # This is OK because the vector is expressed internally as a set of coordinate components in
  # a particular chart, and that representation remains valid as the point moves.

  def __init__(self,point,comp):
    self.point = point
    # During a transition, the following two flags can become different from the ones
    # in the point object:
    self.chart = copy.copy(point.chart)
    self.rot90 = copy.copy(point.rot90)
    if point.chart==SphPoint.KRUSKAL_VW_CHART: self._era=copy.deepcopy(point._era)
    self.comp = comp # array containing four components of the vector, in the current chart
    point.register_vector(self)

  def __str__(self):
    s = self.comp
    return "("+str(s[0])+", "+str(s[1])+", "+str(s[2])+", "+str(s[3])+")"

  def debug_print(self):
    print("  sph_vector.debug_print: ",str(self))

  # This is called automatically from routines in SphPoint for vectors registered as clients of the point.
  # Logically there is the possibility that two things could change at once, e.g., you
  # could change from Schwarzschild to Kruskal, and rot90 or era could also change.
  # The code below can't handle that, and needs to be called multiple times.
  def handle_transition(self):
    if not self.point.transition: return # This should not normally happen.
    chart1 = self.chart
    chart2 = self.point.chart
    if chart1==SphPoint.SCHWARZSCHILD_CHART and chart2==SphPoint.KRUSKAL_VW_CHART:
      v = self.point.v
      w = self.point.w
      t = self.point.t_before_transition
      r = self.point.r_before_transition
      region = schwarzschild.ks_to_region(v,w)
      j = schwarzschild.sch_ks_jacobian(region,r,t)
      dt = self.comp[0] ; dr = self.comp[1]
      dv = j[0][0]*dt+j[0][1]*dr
      dw = j[1][0]*dt+j[1][1]*dr
      self.comp[0]=dv; self.comp[1]=dw
      self.chart = copy.copy(self.point.chart)
      return
    if chart1==SphPoint.KRUSKAL_VW_CHART and chart2==SphPoint.SCHWARZSCHILD_CHART:
      t=self.point.t ; r=self.point.r
      v = self.point.v_before_transition
      w = self.point.w_before_transition
      region = schwarzschild.ks_to_region(v,w)
      j = schwarzschild.ks_sch_jacobian(region,r,t)
      dv = self.comp[0] ; dw = self.comp[1]
      dt = j[0][0]*dv+j[0][1]*dw
      dr = j[1][0]*dv+j[1][1]*dw
      self.comp[0]=dt; self.comp[1]=dr
      self.chart = copy.copy(self.point.chart)
      return
    # Past this point, we're guaranteed that chart1 and chart2 are the same.
    if chart2==SphPoint.KRUSKAL_VW_CHART and self._era!=self.point._era:
      v = self.point.v_before_transition
      w = self.point.w_before_transition
      j = ks_era_jacobian(v,w)
      old_dv = self.comp[0] ; old_dw = self.comp[1]
      dv = j[0][0]*old_dv+j[0][1]*old_dw
      dw = j[0][0]*old_dv+j[0][1]*old_dw
      self.comp[0] = dv ; self.comp[1] = dw
      return
    if self.rot90!=self.point.rot90:
      if (not self.rot90) and self.point.rot90:
        direction = 1.0
      else:
        direction = -1.0
      j = util.jacobian_rot90([self.point.theta_before_transition,self.point.phi_before_transition],direction)
      dtheta = self.comp[2] ; dphi = self.comp[3]
      dtheta2 = j[0][0]*dtheta+j[0][1]*dphi
      dphi2   = j[1][0]*dtheta+j[1][1]*dphi
      self.comp[2]=dtheta2; self.comp[3]=dphi2
    return

  def norm(self):
    g = self.point.metric()
    n = 0.0
    for i in range(4):
      ci = self.comp[i]
      for j in range(4):
        n = n+ci*self.comp[j]*g[i][j]
    return n

  def timelike(self):
    return self.norm()>0.0

