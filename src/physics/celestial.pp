"""
Low-level math functions for coordinates on the celestial sphere and rotation of coordinate systems.
"""

#include "language.h"
#include "util.h"
#include "math.h"
#include "io_util.h"
#include "test.h"
#include "precision.h"

import euclidean

def beta_to_celestial(beta,phi,m):
  """
  An observer is inside the celestial sphere at a certain RA and dec (in radians). They have
  a coordinate system in which beta=0 is the zenith and phi is an azimuthal angle, with phi=0
  being toward the north celestial pole.
  Given a beta and phi, compute the RA and dec of a point below the zenith.
  The matrix m is to have been precomputed using rotation_matrix_observer_to_celestial().
  """
  # Compute the cartesian vector of the point in the observer's frame.
  p = [sin(beta)*cos(phi),sin(beta)*sin(phi),cos(beta)]
  # Rotate to celestial frame:
  p2 = euclidean.apply_matrix(m,p)
  ncp = [0.0,0.0,1.0] # north celestial pole
  cx = [1.0,0.0,0.0] # x axis in celestial coords
  cy = [0.0,1.0,0.0] # x axis in celestial coords
  theta = acos(euclidean.dot(p2,ncp))
  phi = atan2(euclidean.dot(p2,cy),euclidean.dot(p2,cx))
  dec = 0.5*MATH_PI-theta
  ra = copy.copy(phi)#js ra=phi
  if ra<0.0:
    ra = ra + 2.0*MATH_PI
  return [ra,dec]

def celestial_to_beta(ra,dec,m_inv):
  """
  Do the inverse of the transformation done by beta_to_celestial().
  """
  phi = copy.copy(ra)#js phi=ra
  theta = 0.5*MATH_PI-dec
  p = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
  p2 = euclidean.apply_matrix(m_inv,p)
  zenith = [0.0,0.0,1.0]
  ox = [1.0,0.0,0.0] # x axis in obs. coords
  oy = [0.0,1.0,0.0] # x axis in obs. coords
  beta = acos(euclidean.dot(p2,zenith))
  phi = atan2(euclidean.dot(p2,oy),euclidean.dot(p2,ox))
  if phi<0.0:
    phi = phi+2.0*MATH_PI
  return [beta,phi]

def rotation_matrix_observer_to_celestial(ra,dec,direction):
  # See beta_to_celestial() for explanation of what's going on.
  # direction = +-1; if direction is -1, compute the opposite transformation
  # Find the zenith vector, in celestial coords.
  zenith = [cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec)]
  # North celestial pole:
  ncp = [0.0,0.0,1.0]
  # Angle of rotation:
  rot = -direction*acos(euclidean.dot(ncp,zenith))
  # Axis of rotation:
  if rot==0.0:
    axis = ncp # doesn't matter, just avoid division by zero that would otherwise happen in this case
  else:
    axis = euclidean.normalize(euclidean.cross_prod(zenith,ncp))
  return euclidean.rotation_matrix_from_axis_and_angle(rot,axis)

def antipodes_of_ra_and_dec(ra_dec):
  ra,dec = ra_dec
  dec2 = -dec
  ra2 = ra+MATH_PI
  if ra2>2.0*MATH_PI:
    ra2 = ra2 - 2.0*MATH_PI
  return [ra2,dec2]

def lmc_ra_dec():
  """
  Returns the RA and dec of the lesser magellanic cloud, in radians. 
  """
  # https://en.wikipedia.org/wiki/Large_Magellanic_Cloud
  ra = ((5.0+23.0/60.0/3600.0)/24.0)*2*MATH_PI
  dec = (-(69+45/60.0)/360.0)*2*MATH_PI
  return [ra,dec]

def rigel_ra_dec():
  """
  Returns the RA and dec of Rigel, in radians. 
  """
  # 051432.27 -081205.9 +000001.9-000000.600004.2 00.18-0.03B8 0 0.05   2.07, bet Ori, Rigel
  # https://en.wikipedia.org/wiki/Rigel
  ra = ((5.0+14.0/60.0+32.27/3600.0)/24.0)*2*MATH_PI
  dec = (-(8+12/60.0+5.9/3600.)/360.0)*2*MATH_PI
  return [ra,dec]
