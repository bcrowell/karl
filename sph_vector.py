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

  def __init__(self,point,comp):
    self.point = point
    # During a transition, the following two flags can become different from the ones
    # in the point object:
    self.chart = point.chart 
    self.rot90 = point.rot90
    self.comp = comp # array containing four components of the vector, in the current chart

  def __str__(self):
    s = self.comp
    return "("+str(s[0])+", "+str(s[1])+", "+str(s[2])+", "+str(s[3])+")"


  # When updating the point, first set its transition flag to False, then update it, then
  # check whether its transition flag has become true. If so, then call this routine on
  # all vectors that refer to it. There is no reason that this should be called when the
  # chart has not actually changed, so if you do that, it will cause a run-time error
  def handle_transition(self):
    if not self.point.transition: return # This should not normally happen.
    chart1 = self.chart
    chart2 = self.point.chart
    if chart1==SphPoint.SCHWARZSCHILD_CHART and chart2==SphPoint.KRUSKAL_VW_CHART:
      kruskal = self.point.get_raw_coords()
      v=kruskal[0]; w=kruskal[1]
      z = self.point.absolute_schwarzschild() # slow, but we don't care because this happens infrequently
      t = z[0]; r=z[1]
      region = schwarzschild.ks_to_region(v,w)
      j = schwarzschild.sch_ks_jacobian(region,r,t)
      dt = self.comp[0] ; dr = self.comp[1]
      dv = j[0][0]*dt+j[0][1]*dr
      dw = j[1][0]*dt+j[1][1]*dr
      self.comp[0]=dv; self.comp[1]=dw
      return
    if chart1==chart2:
      print("Error, transition from chart ",chart1," to ",chart2," in sph_vector.py, charts are the same")
      raise RuntimeError('')
    print("Error, transition from chart ",chart1," to ",chart2," not yet implemented in sph_vector.py")
    raise RuntimeError('')
    exit(-1)

