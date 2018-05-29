import schwarzschild,util

import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

# Documentation for the math is in the file doc.tex, which can be
# compiled to pdf format by doing a "make doc." (Comments in the code do
# not document the math or the definitions of the variables.)

class SphPoint:
  # This class represents a point in a spherically symmetric spacetime.
  #
  # Instance variables and accessor methods:
  #   spacetime = SCHWARZSCHILD or CHARGED (Reissner-Nordstrom)
  #   chart = SCHWARZSCHILD_CHART or KRUSKAL_VW_CHART
  #   rot90 -- boolean, possible 90-degree rotation to a different chart
  #   v,w = Kruskal null coordinates (defined only if chart is Kruskal)
  #   r,t = Schwarzschild coordinates (defined only if chart if Schwarzschild)
  #   theta,phi = to be used only when we don't care about possible 90-degree rotation
  #   absolute_angles() -- returns [theta,phi] with no arbitrary 90-degree rotation; may be 
  #                        near coordinate singularities at theta=0 and pi; to be used only
  #                        when we don't care about efficiency or coordinate singularities
  #   transition -- boolean, see comments at make_safe()

  # Constants for referring to particular metrics:
  SCHWARZSCHILD = 1
  CHARGED = 2 # Reissner-Nordstrom, i.e., nonzero charge but zero spin
  # Constants referring to different coordinate charts:
  SCHWARZSCHILD_CHART = 1
  KRUSKAL_VW_CHART = 2

  def __init__(self,spacetime,chart,x0,x1,theta,phi):
    self.spacetime = spacetime
    # chart = SCHWARZSCHILD_CHART or KRUSKAL_VW_CHART, the chart in which x0 and x1 are defined
    #   Re Kruskal for Reissner-Nordstrom, see MTW p. 841 (pdf p. 868).
    self.chart = chart
    self.theta = theta
    self.phi = phi
    self.rot90 = False
    if chart==SphPoint.KRUSKAL_VW_CHART:
      self.v = x0
      self.w = x1
      self._era = 0.0
    if chart==SphPoint.SCHWARZSCHILD_CHART:
      self.t = x0
      self.r = x1
    self.make_safe()

  def __str__(self):
    s = self.absolute_schwarzschild()
    return "(t="+str(s[0])+", r="+str(s[1])+", theta="+str(s[2])+", phi="+str(s[3])+")"

  # Add c*v to the point, where v is a vector (must refer to this point), and c is a scalar.
  # After doing this enough times, call make_safe.
  def add(self,c,v):
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      self.v += c*v.comp[0]
      self.w += c*v.comp[1]
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      self.t += c*v.comp[0]
      self.r += c*v.comp[1]
    self.theta += c*v.comp[2]
    self.phi   += c*v.comp[3]

  # This is for use when calculating things like christoffel symbols, where we don't
  # care about rot90 or era because they don't affect the results.
  def get_raw_coords(self):
    coords = [0.0,0.0,0.0,0.0]
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      coords[0] = self.v
      coords[1] = self.w
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      coords[0] = self.t
      coords[1] = self.r
    coords[2] = self.theta
    coords[3] = self.phi
    return coords

  def set_raw_coords(self,coords):
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      self.v = coords[0]
      self.w = coords[1]
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      self.t = coords[0]
      self.r = coords[1]
    self.theta = coords[2]
    self.phi = coords[3]

  # Public routine that does manipulations on the internal representation in order to avoid
  # the following three issues:
  #   - coordinate singularities at theta=0 and pi
  #   - numerical misbehavior of Kruskal coordinates for large r and t
  #   - coordinate singularity in the Schwarzschild chart at the event horizon
  # External code should call this method often enough so that it happens
  # before:
  #   - the angles can change by ~0.1
  #   - Kruskal V,W can change by ~1
  #   - motion can occur from large r>>r_Sch to the horizon
  # If there are tangent vectors expressed using this point, then it may be necessary to
  # reexpress them in a different chart after calling this routine. To detect this, set
  # the point's .transition flag to False before calling, and then check it afterward.
  def make_safe(self):
    self._rotate_to_safety()
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      self._era_in_range()
      rho = -self.v*self.w
      if rho>2000.0: self.to_schwarzschild()
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      if self.r<3.0: self.to_kruskal()

  def to_schwarzschild(self):
    if self.chart==SphPoint.SCHWARZSCHILD_CHART: return
    tr = schwarzschild.ks_to_sch(self.v,self.w)
    self.t = tr[0]
    self.r = tr[1]
    self.chart = SphPoint.SCHWARZSCHILD_CHART
    self.transition = True

  def to_kruskal(self):
    if self.chart==SphPoint.KRUSKAL_VW_CHART: return
    vw = schwarzschild.sch_to_ks(self.t,self.r,schwarzschild.sch_to_sigma(self.r))
    self.v = vw[0]
    self.w = vw[1]
    self._era = 0.0
    self.chart = SphPoint.KRUSKAL_VW_CHART
    self._era_in_range()
    self.transition = True

  def absolute_angles(self):
    angles = [self.theta,self.phi]
    if self.rot90: angles = util.rotate_unit_sphere(angles,-1.0)
    return angles

  # Return an array containing Kruskal null coordinates [v,w,theta,phi], in absolute
  # form for output. This routine is not meant to be used for internal manipulations.
  # E.g., when calculating a Christoffel coefficient, just use .v, etc., to exploit
  # the rotational and time symmetries, and to avoid overflows. Note that if r or t is large,
  # this may cause a floating-point exception.
  def absolute_kruskal(self):
    angles = [self.theta,self.phi]
    if self.rot90: angles = util.rotate_unit_sphere(angles,-1.0)
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      vw = schwarzschild.ks_to_zero_era(self._era,self.v,self.w)
    else:
      vw = schwarzschild.sch_to_ks(self.t,self.r,schwarzschild.sch_to_sigma(self.r))
    return [vw[0],vw[1],angles[0],angles[1]]

  # Similar to absolute_kruskal. Slow, intended to be used for output.
  def absolute_schwarzschild(self):
    angles = [self.theta,self.phi]
    if self.rot90: angles = util.rotate_unit_sphere(angles,-1.0)
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      tr = schwarzschild.ks_to_sch(self.v,self.w)
      t = tr[0]+self._era
      r = tr[1]
    else:
      t = self.t
      r = self.r
    return [t,r,angles[0],angles[1]]

  # Do not call this routine unless the chart is already known to be Kruskal.
  def _era_in_range(self):
    ks = schwarzschild.ks_era_in_range([self._era,self.v,self.w])
    self._era = ks[0]
    self.v = ks[1]
    self.w = ks[2]
    self.transition = True

  # The colatitude theta is supposed to be in [0,pi]. If it's slightly out
  # of that range, bump it in. This routine assumes that if we're out of
  # range, it's due to continuous motion, so it doesn't try very hard if
  # we're way out of range. Similar manipulations for phi.
  def _canonicalizetheta(self):
    if self.theta<0.0: self.theta = -self.theta
    if self.theta>pi:
      self.theta = 2.0*pi-self.theta
      self.phi = self.phi + 2.0*pi
    if self.phi>2.0*pi: self.phi -= 2.0*pi
    if self.phi<0.0: self.phi += 2.0*pi

  # If necessary, rotate into a different coordinate chart in order to avoid the coordinate
  # singularities at theta=0 and pi.
  def _rotate_to_safety(self):
    self._canonicalizetheta()
    if self.theta<0.3 or self.theta>2.841: # pi-0.3
      if self.rot90:
        direction = -1.0
      else:
        direction = 1.0
      self.rot90 = not self.rot90
      angles = util.rotate_unit_sphere([self.theta,self.phi],direction)
      self.theta = angles[0]
      self.phi = angles[1]

