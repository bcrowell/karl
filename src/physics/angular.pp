"""
Helper routines to compute things for the representation of the
unit sphere as a sphere embedded in three-dimensional cartesian
space (i,j,k).
"""

#include "language.h"
#include "util.h"
#include "math.h"

def renormalize(x):
  """
  Given a point represented by coordinates of the form (...,...,i,j,k),
  return a copy with i, j, and k renormalized to lie on the unit sphere.
  """
  n = CLONE_ARRAY_OF_FLOATS(x)
  xi = sqrt(n[2]*n[2]+n[3]*n[3]+n[4]*n[4])
  n[2] = n[2]/xi
  n[3] = n[3]/xi
  n[4] = n[4]/xi
  return n

def make_tangent(x,v0):
  """
  Given a point x and a velocity vector v in its tangent space, return a copy of
  v in which the component perpendicular to the i-j-k sphere has been removed.

  For accurate results, x should already be accurately normalized. Since only the
  angular parts of x and v0 are used, it doesn't matter whether they are expressed
  in Schwarzschild or Kruskal, as long as they're 5-dimensional.
  """
  v = CLONE_ARRAY_OF_FLOATS(v0)
  dot = x[2]*v[2]+x[3]*v[3]+x[4]*v[4]
  v[2] -= dot*x[2]
  v[3] -= dot*x[3]
  v[4] -= dot*x[4]
  return v

def theta_phi_to_ijk(theta,phi):
  sin_theta = sin(theta)
  return [sin_theta*cos(phi),sin_theta*sin(phi),cos(theta)]

def add_centripetal(ch,p):
  """
  Add centripetal terms to the Christoffel symbols for a 5-dimensional coordinate system.
  Modifies ch in place by adding the centripetal parts.
  """
  i=p[2]
  j=p[3]
  k=p[4]
  xi2 = i*i+j*j+k*k; # should normally be very close to 1
  for m in range(2,5): # upper index
    z = p[m]
    for n in range(2,5): # lower indices
      ch[n][n][m] += z/xi2
