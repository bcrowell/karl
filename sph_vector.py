import copy

import schwarzschild,util

import numpy as np
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
    self.chart = copy.deepcopy(point.chart)
    self.rot90 = copy.deepcopy(point.rot90)
    if point.chart==SphPoint.KRUSKAL_VW_CHART: self._era=copy.deepcopy(point._era)
    self.comp = comp # array containing four components of the vector, in the current chart

  def __str__(self):
    s = self.comp
    return "("+str(s[0])+", "+str(s[1])+", "+str(s[2])+", "+str(s[3])+")"


  # When updating the point, first set its transition flag to False, then update it, then
  # check whether its transition flag has become true. If so, then call this routine on
  # all vectors that refer to it. There is no reason that this should be called when nothing
  # has actually changed, so if you do that, it is not guaranteed to work and
  # may cause a run-time error. (In any case, it will be inefficient to do so.)
  # Logically there is the possibility that two things could change at once, e.g., you
  # could change from Schwarzschild to Kruskal, and rot90 or era could also change.
  # In fact the only combination that is allowed to occur is if both rot90 and era change.
  # FIXME -- this is not OK, change to Kruskal can also trigger era.
  # We need the info about the old chart and coordinates, which are preserved in fields
  # such as t_before_transition.
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
      return # don't fall through, because nothing else is allowed to change if the chart changed
    if chart1==SphPoint.KRUSKAL_VW_CHART and chart2==SphPoint.SCHWARZSCHILD_CHART:
      t=self.point.t ; r=self.point.r
      v = self.point.v_before_transition
      w = self.point.w_before_transition
      region = ks_to_region(v,w)
      j = schwarzschild.ks_sch_jacobian(region,r,t)
      dv = self.comp[0] ; dw = self.comp[1]
      dt = j[0][0]*dv+j[0][1]*dw
      dr = j[1][0]*dv+j[1][1]*dw
      self.comp[0]=dt; self.comp[1]=dr
      return # don't fall through, because nothing else is allowed to change if the chart changed
    # Past this point, we're guaranteed that chart1 and chart2 are the same.
    if chart2==SphPoint.KRUSKAL_VW_CHART and self._era!=self.point._era:
      v = self.point.v_before_transition
      w = self.point.w_before_transition
      j = ks_era_jacobian(v,w)
      old_dv = self.comp[0] ; old_dw = self.comp[1]
      dv = j[0][0]*old_dv+j[0][1]*old_dw
      dw = j[0][0]*old_dv+j[0][1]*old_dw
      self.comp[0] = dv ; self.comp[1] = dw
    if self.rot90!=self.point.rot90:
      if (not self.rot90) and self.point.rot90:
        direction = 1.0
      else:
        direction = -1.0
      old_angles = util.rotate_unit_sphere([self.comp[2],self.comp[3]],-direction) # find what angles were before transition
      theta = old_angles[0]; phi = old_angles[1]
      j = util.jacobian_rot90([theta,phi],direction)
      dtheta = self.comp[2] ; dphi = self.comp[1]
      dtheta2 = j[0][0]*dtheta+j[0][1]*dphi
      dphi2   = j[1][0]*dtheta+j[1][1]*dphi
      self.comp[2]=dtheta2; self.comp[3]=dphi2





