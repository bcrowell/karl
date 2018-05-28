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
  # all vectors that refer to it.
  def check_for_transition(self):
    if not self.point.transition: return # This should not normally happen.
    print("Error, transition code not yet implemented in sph_vector.py")
    exit(-1)

