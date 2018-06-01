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

# The jacobian matrix of the transformation described by rotate_unit_sphere().
# This is nearly raw output from the following Maxima code:
#   direction:1;
#   s:sin(theta);
#   theta2:acos(-direction*s*sin(phi));
#   phi2:atan2(direction*cos(theta),s*cos(phi));
#   j:jacobian([theta2,phi2],[theta,phi])$
#   print("1 1 ",string(j[1,1]));
#   ...etc
# This could be optimized a lot more, but it's called infrequently so that's not really worthwhile.
# The order of the subscripts is such that, e.g., j[1][0] is dphi2/dtheta.
def jacobian_rot90(angles,direction):
  theta = angles[0]
  phi = angles[1]
  j = [[0 for i in range(2)] for j in range(2)]
  a = sqrt(1-sin(phi)**2*sin(theta)**2)
  if a==0.0:
    raise RuntimeError("error in util.jacobian_rot90, a=0, phi="+str(phi)+", theta="+str(theta))
    # This shouldn't happen, because we transition before we enter the polar "cap."
  if direction>0.0:
    j[0][0] = (sin(phi)*cos(theta))/a
    j[0][1] = (cos(phi)*sin(theta))/a 
    j[1][0] = (-(cos(phi)*sin(theta)**2)/(cos(theta)**2+cos(phi)**2*sin(theta)**2))-(cos(phi)*cos(theta)**2)/(cos(theta)**2+cos(phi)**2*sin(theta)**2) 
    j[1][1] = (sin(phi)*cos(theta)*sin(theta))/(cos(theta)**2+cos(phi)**2*sin(theta)**2) 
  else:
    j[0][0] = -(sin(phi)*cos(theta))/a 
    j[0][1] = -(cos(phi)*sin(theta))/a 
    j[1][0] = (cos(phi)*sin(theta)**2)/(cos(theta)**2+cos(phi)**2*sin(theta)**2)+(cos(phi)*cos(theta)**2)/(cos(theta)**2+cos(phi)**2*sin(theta)**2) 
    j[1][1] = -(sin(phi)*cos(theta)*sin(theta))/(cos(theta)**2+cos(phi)**2*sin(theta)**2) 
  return j

