"""
Low-level routines to compute transformations between Schwarzschild and
arcsinh-Kruskal coordinates, and the jacobians of those transformations.

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.) 
"""

#include "language.h"
#include "util.h"
#include "math.h"
#include "precision.h"

import math_util
import kruskal

def schwarzschild_to_kruskal(t,r):
  """
  Transforms Schwarzschild (t,r) coordinates to arcsinh-Kruskal (a,b). Assumes region I or II.
  """
  p = sqrt(abs(r-1))
  if r>1.0:
    sign_b = -1.0
  else:
    sign_b = 1.0
  return [schwarzschild_to_kruskal_helper(p,r+t),sign_b*schwarzschild_to_kruskal_helper(p,r-t)]

def schwarzschild_to_kruskal_helper(p,q):
  # compute y=asinh(pe^(q/2))
  return math_util.asinh_of_exp(log(p)+q/2)

def kruskal_to_schwarzschild(a,b):
  """
  Transforms arcsinh-Kruskal (a,b) coordinates to Schwarzschild (t,r).

  This is really just a wrapper for kruskal.aux().
  """
  t,r,mu = kruskal.aux(a,b)
  return [t,r]

#-----------------------------------------------------------------------------------------------------

def jacobian_schwarzschild_to_kruskal(t,r):
  """
  Returns the matrix of partial derivatives of (a,b) with respect to (t,r), given t and r.
  The first index is 0 for a, 1 for b. The second index is 0 for t, 1 for r.

  For non-horizon points, assumes region I or II.
  For points on the horizon, the result contains some infinite matrix elements.
  """
  jacobian = EMPTY2DIM(2)
  log_r_minus_1 = log(abs(r-1))
  xpi2 = math_util.safe_exp(-log_r_minus_1-(r+t)) # x_+^{-2}
  xmi2 = math_util.safe_exp(-log_r_minus_1-(r-t)) # x_-^{-2}
  if r>1.0:
    s = 1.0 # region I
  else:
    s = -1.0 # region II
  jacobian[0][0] =   0.5/sqrt(1+xpi2) # da/dt
  jacobian[1][0] = s*0.5/sqrt(1+xmi2) # db/dt
  if r!=1.0:
    # not on the horizon
    q = 1.0/(1.0-1.0/r)
    jacobian[0][1] =   q*jacobian[0][0] # da/dr
    jacobian[1][1] =  -q*jacobian[1][0] # db/dr
  else:
    # on the horizon
    jacobian[0][1] = float("inf")
    jacobian[1][1] = float("inf")
  return jacobian  

def jacobian_kruskal_to_schwarzschild(t,r):
  """
  Returns the matrix of partial derivatives of (t,r) with respect to (a,b), given t and r.
  The first index is 0 for t, 1 for r. The second index is 0 for a, 1 for b.

  Slightly inefficient because we invert the 2x2 matrix,
  but not likely to affect performance because not called often.
  This can throw an error if called for a point on the horizon (r=1), where
  the Schwarzschild coordinates misbehave, and in general it's probably
  not going to be numerically accurate to use this near the horizon.
  """
  if r==1.0:
    THROW('r=1 in jacobian_kruskal_to_schwarzschild')
  j1 = jacobian_schwarzschild_to_kruskal(t,r)
  a = j1[0][0]
  d = j1[1][1]
  b = j1[1][0]
  c = j1[0][1]
  det = a*d-b*c # should be nonzero because we checked above for r==1
  j2 = EMPTY2DIM(2) # allocate a new matrix, since a-d are just pointers into j1
  j2[1][1] = a/det
  j2[0][0] = d/det
  j2[1][0] = -b/det
  j2[0][1] = -c/det
  return j2

#-----------------------------------------------------------------------------------------------------

def schwarzschild_to_kruskal_small(t,r):
  """
  Compute Kruskal (a,b) coordinates from Schwarzschild (t,r). Overflows if t and r are not small.
  Returns a result in region I or II. For testing purposes only.
  """
  t2 = 0.5*t
  if r>1.0:
    sc = sinh(t2)
    cs = cosh(t2)
  else:
    sc = cosh(t2)
    cs = sinh(t2)
  # My formulation based on MTW p. 827:
  h = sqrt(abs(r-1.0))*exp(r/2.0)
  ks_t = h*sc
  ks_x = h*cs
  v = ks_t+ks_x
  w = ks_t-ks_x
  return [arcsinh(v),arcsinh(w)]
