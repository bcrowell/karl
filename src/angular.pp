"""
Helper routines to compute things for the representation of the
unit sphere as a sphere embedded in three-dimensional cartesian
space (i,j,k).
"""

#include "util.h"
#include "math.h"

def renormalize(x):
  """
  Given a point represented by coordinates of the form (...,...,i,j,k),
  return a copy with i, j, and k renormalized to lie on the unit sphere.
  """
  n = copy.copy(x)
  xi = sqrt(n[2]*n[2]+n[3]*n[3]+n[4]*n[4])
  n[2] = n[2]/xi
  n[3] = n[3]/xi
  n[4] = n[4]/xi
  return n

def make_tangent(x,v0):
  """
  Given a point x and a velocity vector v in its tangent space, return a copy of
  v in which the component perpendicular to the i-j-k sphere has been removed.

  For accurate results, x should already be accurately normalized.
  """
  v = copy.copy(v0)
  dot = x[2]*v[2]+x[3]*v[3]+x[4]*v[4]
  v[2] -= dot*x[2]
  v[3] -= dot*x[3]
  v[4] -= dot*x[4]
  return v

