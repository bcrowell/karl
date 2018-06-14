"""
Low-level routines to compute transformations between Schwarzschild and
arcsinh-Kruskal coordinates, and the jacobians of those transformations.

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.) 
"""

#include "util.h"
#include "math.h"
#include "precision.h"

import math_util
from math_util import safe_exp

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

#-----------------------------------------------------------------------------------------------------

def jacobian_schwarzschild_to_kruskal(t,r):
  """
  Returns the matrix of partial derivatives of (a,b) with respect to (t,r). 
  The first index is 0 for a, 1 for b. The second index is 0 for t, 1 for r.

  For non-horizon points, assumes region I or II.
  For points on the horizon, the result contains some infinite matrix elements.
  """
  jacobian = [[0 for i in range(2)] for j in range(2)]
  log_r_minus_1 = log(abs(r-1))
  xpi2 = safe_exp(-log_r_minus_1-(r+t)) # x_+^{-2}
  xmi2 = safe_exp(-log_r_minus_1-(r-t)) # x_-^{-2}
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
