from scipy.special import lambertw

import numpy as np
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi
from scipy import sign

###################################################################

def lambert_w(x):
  return lambertw(x).real

###################################################################

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
