import copy

import schwarzschild,util,io_util

import numpy
numpy.seterr(all='raise')
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
  #   chart_before_transition -- what the chart was before the most recent transition
  #   get_christoffel() -- returns the christoffel symbols in the current coordinate chart
  #   sch_copy --- when using Schwarzschild coordinates, this keeps track of the region
  #                of the maximal extension of the Schwarzschild spacetime; if this boolean flag
  #                is True, then we're in region 3 or 4

  # Constants for referring to particular metrics:
  SCHWARZSCHILD = 1
  CHARGED = 2 # Reissner-Nordstrom, i.e., nonzero charge but zero spin
  # Constants referring to different coordinate charts:
  SCHWARZSCHILD_CHART = 1
  KRUSKAL_VW_CHART = 2
  # Constants that control when we transition from one chart to another:
  TRANSITION_MAX_RHO = 2000.0 # max value of rho before we transition from Kruskal to Sch.
  TRANSITION_MIN_R = 3.0      # min value of r before we transition from Sch. to Kruskal
  TRANSITION_THETA = 0.3      # how close to poles before doing a rot90

  # After creating a new point with this code, you may immediately want to call its make_safe()
  # method to get into a more appropriate coordinate chart. However, often you will create
  # vectors in the tangent space associate with the point as well. E.g., if simulating a geodesic,
  # you will first create an initial position x, and then create an initial tangent vector v.
  # You want to create the tangent vector in the expected chart, so you would typically do that,
  # and *then* if necessary transition both x and v to the new chart. This use case is the reason why
  # make_safe() is not automatically called when a new point is created. Because make_safe() is
  # not called automatically, the .transition flag is guaranteed never to be set after creating a new point.
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
      self.sch_copy = False
    self.client_vectors = []
    self.debug_transitions = False

  def __str__(self):
    describe_angles = "theta="+("%5.3e" % self.theta)+", phi="+("%5.3e" % self.phi)
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      return "(V="+("%5.3e" % self.v)+", W="+("%5.3e" % self.w)+", "+describe_angles+")"
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      return "(t="+("%5.3e" % self.t)+", r="+("%5.3e" % self.r)+", "+describe_angles+")"
    return "error stringifying a SphPoint"

  def debug_print(self,header):
    print("==== ",header," ===")
    print("sph_point.debug_print: ")
    print("  chart = ",self.chart)
    print("  coordinates: ",str(self))
    print("  Schwarzschild coordinates: ",io_util.vector_to_str(self.absolute_schwarzschild()))
    print("  number of registered vectors = ",len(self.client_vectors))
    self.debug_print_vectors()

  def debug_print_vectors(self):
    print("debug_print_vectors:")
    for v in self.client_vectors:
      v.debug_print()

  def register_vector(self,v):
    self.client_vectors.append(v)

  def deregister_vector(self,victim):
    self.client_vectors = [v for v in self.client_vectors if v != victim]

  def metric(self):
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      return schwarzschild.sch_ks_metric(self.v,self.w,self.theta)
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      return schwarzschild.sch_ks_metric(self.r,self.theta)
    raise RuntimeError('unrecognized chart')

  def region(self):
    if self.chart==SphPoint.KRUSKAL_VW_CHART: return ks_to_region(self.v,self.w)
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      if r<1.0:
        if self.sch_copy:
          return 4
        else:
          return 2
      else:
        if self.sch_copy:
          return 3
        else:
          return 1
    raise RuntimeError('unknown chart in sph_point.py, region()')

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
    self.transition = False
    self._rotate_to_safety()
    if self.chart==SphPoint.KRUSKAL_VW_CHART:
      self._era_in_range()
      rho = -self.v*self.w
      if rho>SphPoint.TRANSITION_MAX_RHO: self.to_schwarzschild()
    if self.chart==SphPoint.SCHWARZSCHILD_CHART:
      # If changing the parameter 3.0, then change the test in test_geodesic_rk_free_fall_from_rest so we
      # still have coverage.
      if self.r<SphPoint.TRANSITION_MIN_R: self.to_kruskal()

  def transition_vectors_if_necessary(self):
    if not self.transition: return
    for v in self.client_vectors:
      v.handle_transition()
    self.transition = False

  # This will have the side-effect of setting the .transition flag.
  def force_chart(self,desired_chart):
    if self.chart==desired_chart: return
    if desired_chart==SphPoint.SCHWARZSCHILD_CHART:
      self.to_schwarzschild()
      return
    if desired_chart==SphPoint.KRUSKAL_VW_CHART:
      self.to_kruskal()
      return
    raise RuntimeError("chart "+str(desired_chart)+" not implemented in sph_point.force_chart")

  def to_schwarzschild(self):
    if self.chart==SphPoint.SCHWARZSCHILD_CHART: return
    self.chart_before_transition = self.chart
    self.v_before_transition = copy.deepcopy(self.v)
    self.w_before_transition = copy.deepcopy(self.w)
    tr = schwarzschild.ks_to_sch(self.v,self.w)
    self.t = tr[0]
    self.r = tr[1]
    self.chart = SphPoint.SCHWARZSCHILD_CHART
    self.sch_copy = (self.v<0.0)
    self.transition = True
    self.transition_vectors_if_necessary()
    if self.debug_transitions: self.debug_print("did a change of chart from Kruskal to Sch.")

  def to_kruskal(self):
    if self.chart==SphPoint.KRUSKAL_VW_CHART: return
    if self.debug_transitions: self.debug_print("doing a change of chart from Sch. to Kruskal")
    self.chart_before_transition = self.chart
    self.t_before_transition = copy.deepcopy(self.t)
    self.r_before_transition = copy.deepcopy(self.r)
    vw = schwarzschild.sch_to_ks(self.t,self.r,schwarzschild.sch_to_sigma(self.r))
    self.v = vw[0]
    self.w = vw[1]
    self._era = 0.0
    self.chart = SphPoint.KRUSKAL_VW_CHART
    self._era_in_range()
    self.transition = True
    self.transition_vectors_if_necessary()
    if self.debug_transitions: self.debug_print("did a change of chart from Sch. to Kruskal")

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
    self.v_before_transition = copy.deepcopy(self.v)
    self.w_before_transition = copy.deepcopy(self.w)
    z = schwarzschild.ks_era_in_range([self._era,self.v,self.w])
    if not z[0]: return # no change was needed
    ks = z[1]
    self.chart_before_transition = self.chart
    self._era = ks[0]
    self.v = ks[1]
    self.w = ks[2]
    self.transition = True
    self.transition_vectors_if_necessary()
    if self.debug_transitions: self.debug_print("did a change of era")

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
    if self.theta>SphPoint.TRANSITION_THETA and self.theta<pi-SphPoint.TRANSITION_THETA: return # [0.3,pi-0.3]
    self.theta_before_transition = copy.copy(self.theta)
    self.phi_before_transition = copy.copy(self.phi)
    if self.debug_transitions: self.debug_print("about to do a 90-degree rotation")
    if self.rot90:
      direction = -1.0
    else:
      direction = 1.0
    self.rot90 = not self.rot90
    angles = util.rotate_unit_sphere([self.theta,self.phi],direction)
    self.theta = angles[0]
    self.phi = angles[1]
    if self.debug_transitions: self.debug_print("rotated point")
    self.transition = True
    self.transition_vectors_if_necessary()
    if self.debug_transitions: self.debug_print("did a 90-degree rotation")

  # Get the Christoffel symbols in the current coordinate chart.
  # In principle the coords argument is superfluous, because we have the data internally
  # as part of the SphPoint object. However, we are normally going to be calling this during
  # Runge-Kutta integration, where we don't want the overhead of creating lots of modified
  # copies of the original point.
  def get_christoffel(self,coords):
    theta = coords[2] # we don't care about rot90 in this context
    if self.spacetime==SphPoint.SCHWARZSCHILD:
      if self.chart==SphPoint.SCHWARZSCHILD_CHART:
        return schwarzschild.sch_christoffel_sch(coords[0],coords[1],sin(theta),cos(theta))
      if self.chart==SphPoint.KRUSKAL_VW_CHART:
        z = schwarzschild.sch_aux_ks(coords[0],coords[1])
        return schwarzschild.sch_christoffel_ks(coords[0],coords[1],sin(theta),cos(theta),z[1],z[2])
    print("error in sph_point.get_christoffel, unimplemented spacetime ",self.spacetime," or chart ",self.chart)
    exit(-1)

