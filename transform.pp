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

def schwarzschild_to_kruskal(t,r):
  p = sqrt(abs(r-1))
  if r>1.0:
    sign_b = -1.0
  else:
    sign_b = 1.0
  return [from_schwarzschild_helper(p,r+t),sign_b*from_schwarzschild_helper(p,r-t)]

def from_schwarzschild_helper(p,q):
  # compute y=asinh(pe^(q/2))
  return math_util.asinh_of_exp(log(p)+q/2)

def from_schwarzschild_small(t,r):
  """
  Compute Kruskal (a,b) coordinates from Schwarzschild (t,r). Overflows if t and r are not small.
  Returns a result in region I or II.
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
