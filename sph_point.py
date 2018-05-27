import schwarzschild

import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

class SphPoint:
  # This class represents a point in a spherically symmetric spacetime.

  # Constants for referring to particular metrics:
  SCHWARZSCHILD = 1
  CHARGED = 2 # Reissner-Nordstrom, i.e., nonzero charge but zero spin
  # Constants referring to different coordinate charts:
  SCHWARZSCHILD_CHART = 1
  KRUSKAL_VW = 2

  def __init__(self,spacetime,chart,x0,x1,theta,phi):
    self.spacetime = spacetime # SCHWARZSCHILD or CHARGED
    # chart = SCHWARZSCHILD_CHART or KRUSKAL_VW, the chart in which x0 and x1 are defined
    #   Re Kruskal for Reissner-Nordstrom, see MTW p. 841 (pdf p. 868).
    self.theta = theta
    self.phi = phi
    ok = False
    if chart==SphPoint.KRUSKAL_VW:
      self.v = x0
      self.w = x1
      ok = True
    if chart==SphPoint.SCHWARZSCHILD_CHART and spacetime==SphPoint.SCHWARZSCHILD:
      t = x0
      r = x1
      vw = schwarzschild.sch_to_ks(t,r,schwarzschild.sch_to_sigma(r))
      self.v = vw[0]
      self.w = vw[1]
      ok = True
    if not ok: die("Error creating SphPoint object, Kruskal-like coordinates not implemented for Reissner-Nordstrom")
    self._era = 0.0
    self._theta = theta
    self._phi = phi
    self._rot90 = False
    self.make_safe()

  # Public routine that does manipulations on the internal representation in order to avoid
  # coordinate singularities at the poles and numerical misbehavior of Kruskal coordinates
  # for large times. External code should call this code often enough so that it happens
  # before the angles can change by ~0.1 or Kruskal V,W can change by ~1.
  def make_safe(self):
    self._era_in_range()
    self._rotate_to_safety()

  # Return an array containing Kruskal null coordinates [v,w,theta,phi], in absolute
  # form for output. This routine is not meant to be used for internal manipulations.
  # E.g., when calculating a Christoffel coefficients, just use .v, etc., to exploit
  # the rotational and time symmetries, and to avoid overflows. Note that if r or t is large,
  # this may cause a floating-point exception.
  def absolute_kruskal(self):
    angles = [self._theta,self._phi]
    if self._rot90: angles = rotate_unit_sphere(angles,-1.0)
    vw = schwarzschild.ks_to_zero_era(self._era,self.v,self.w)
    return [vw[0],vw[1],angles[0],angles[1]]

  def _era_in_range(self):
    ks = schwarzschild.ks_era_in_range([self._era,self.v,self.w])
    self._era = ks[0]
    self.v = ks[1]
    self.w = ks[2]

  # The colatitude theta is supposed to be in [0,pi]. If it's slightly out
  # of that range, bump it in. This routine assumes that if we're out of
  # range, it's due to continuous motion, so it doesn't try very hard if
  # we're way out of range. Similar manipulations for phi.
  def _canonicalize_theta(self):
    if self._theta<0.0: self._theta = -self._theta
    if self._theta>pi:
      self._theta = 2.0*pi-self._theta
      self._phi = self._phi + 2.0*pi
    if self._phi>2.0*pi: self._phi -= 2.0*pi
    if self._phi<0.0: self._phi += 2.0*pi

  # If necessary, rotate into a different coordinate chart in order to avoid the coordinate
  # singularities at theta=0 and pi.
  def _rotate_to_safety(self):
    self._canonicalize_theta()
    if self._theta<0.3 or self.theta>2.841: # pi-0.3
      if self._rot90:
        direction = -1.0
      else:
        direction = 1.0
      self._rot90 = not self._rot90
      angles = rotate_unit_sphere([self._theta,self._phi],direction)
      self._theta = angles[0]
      self._phi = angles[1]


  def die(self,message):
    print(message)
    exit(-1)
