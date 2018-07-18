"""
Routines to compute transformations of points among Schwarzschild,
arcsinh-Kruskal, and Keplerian coordinates, the jacobians of those transformations,
and transformations of vectors.

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.) 
"""

#include "language.h"
#include "util.h"
#include "math.h"
#include "precision.h"
#include "spacetimes.h"

import math_util
import kruskal,keplerian,schwarzschild,io_util

def chart_info(spacetime_and_chart):
  """
  Input is spacetime|chart. Returns [ok,ndim,christoffel_function,name].

  See README for a list of things to do when adding a new chart.
  """
  recognized = FALSE
  if (spacetime_and_chart)==(SP_SCH|CH_SCH):
    return [TRUE,5,schwarzschild.christoffel,"SCH"]
  if (spacetime_and_chart)==(SP_SCH|CH_AKS):
    return [TRUE,5,kruskal.christoffel,"AKS"]
  if (spacetime_and_chart)==(SP_SCH|CH_KEP):
    return [TRUE,5,keplerian.christoffel,"KEP"]
  return [FALSE,None,None,'']

def transform_point(x,spacetime,chart,pars,chart2):
  """
  Transforms a point x from chart to chart2. Return value is an array, which is automatically cloned.
  It's legal to have chart==chart2. If the input point is not in the domain of the relevant functions,
  returns an array whose elements are all NaN.
  """
  if chart==chart2:
    return CLONE_ARRAY_OF_FLOATS(x)
  if spacetime!=SP_SCH:
    THROW_ARRAY((["unrecognized spacetime=",spacetime]))
  ok,ndim,christoffel_function,name = chart_info(spacetime|chart)
  if not ok:
    THROW_ARRAY((["unrecognized chart, spacetime=",spacetime,"chart=",chart]))
  ok,ndim2,christoffel_function,name = chart_info(spacetime|chart2)
  if not ok:
    THROW_ARRAY((["unrecognized chart, spacetime=",spacetime,"chart=",chart2]))
  x2=EMPTY1DIM(ndim2)
  if (chart==CH_SCH or chart==CH_KEP) and x[1]<0.0:
    for i in range(ndim2):
      x2[i] = NAN
    return CLONE_ARRAY_OF_FLOATS(x2) # cloning not really needed in this case, but be consistent
  done = FALSE
  # The following assumes, as is true for CH_SCH, CH_AKS, and CH_KEP, that coords 2,3,4 are (i,j,k).
  x2[2]=x[2]
  x2[3]=x[3]
  x2[4]=x[4]
  # Transform the first two coordinates:
  if chart==CH_SCH and chart2==CH_AKS:
    a,b = schwarzschild_to_kruskal(x[0],x[1])
    x2[0] = a
    x2[1] = b
    done = TRUE
  if not done and chart==CH_AKS and chart2==CH_SCH:
    t,r = kruskal_to_schwarzschild(x[0],x[1])
    x2[0] = t
    x2[1] = r
    done = TRUE
  if not done and chart==CH_SCH and chart2==CH_KEP:
    x2[0] = x[0]
    x2[1] = x[1]**(3.0/2.0)
    done = TRUE
  if not done and chart==CH_KEP and chart2==CH_SCH:
    x2[0] = x[0]
    x2[1] = x[1]**(2.0/3.0)
    done = TRUE
  if not done and chart==CH_AKS and chart2==CH_KEP:
    t,r = kruskal_to_schwarzschild(x[0],x[1])
    x2[0] = t
    x2[1] = r**(3.0/2.0)
    done = TRUE
  if not done and chart==CH_KEP and chart2==CH_AKS:
    t = x[0]
    r = x[1]**(2.0/3.0)
    a,b = schwarzschild_to_kruskal(t,r)
    x2[0] = a
    x2[1] = b
    done = TRUE
  if done:
    return CLONE_ARRAY_OF_FLOATS(x2) # for python, divorce the new vector from entanglement with components of old
  else:
    THROW_ARRAY((["don't know how to transform, spacetime=",spacetime,", chart=",chart,", chart2=",chart2]))

def transform_vector(v,x,spacetime,chart,pars,chart2):
  """
  Transforms a vector v from chart to chart2 in the tangent space at x (x being described in chart, not chart2).
  Return value is an array, which is
  automatically cloned. It's legal to have chart==chart2.
  """
  if chart==chart2:
    return CLONE_ARRAY_OF_FLOATS(v)
  if spacetime!=SP_SCH:
    THROW_ARRAY((["unrecognized spacetime=",spacetime]))
  ok,ndim,christoffel_function,name = chart_info(spacetime|chart)
  if not ok:
    THROW_ARRAY((["unrecognized chart, spacetime=",spacetime,"chart=",chart]))
  ok,ndim2,christoffel_function,name = chart_info(spacetime|chart2)
  if not ok:
    THROW_ARRAY((["unrecognized chart, spacetime=",spacetime,"chart=",chart2]))
  v2=EMPTY1DIM(ndim2)
  # The following assumes, as is true for CH_SCH, CH_AKS, and CH_KEP, that coords 2,3,4 are (i,j,k).
  v2[2]=v[2]
  v2[3]=v[3]
  v2[4]=v[4]
  found_jac = FALSE
  # Find the Jacobian.
  if chart==CH_SCH and chart2==CH_AKS:
    t,r = [x[0],x[1]]
    jac = jacobian_schwarzschild_to_kruskal(t,r)
    found_jac = TRUE
  if chart==CH_AKS and chart2==CH_SCH:
    t,r = kruskal_to_schwarzschild(x[0],x[1])
    jac = jacobian_kruskal_to_schwarzschild(t,r)
    found_jac = TRUE
  if chart==CH_SCH and chart2==CH_KEP:
    t,r = [x[0],x[1]]
    jac = EMPTY2DIM(2)
    jac[0][0] = 1.0
    jac[1][0] = 0.0
    jac[0][1] = 0.0
    jac[1][1] = (1.5)*sqrt(r)
    found_jac = TRUE
  if chart==CH_KEP and chart2==CH_SCH:
    t,u = [x[0],x[1]]
    jac = EMPTY2DIM(2)
    jac[0][0] = 1.0
    jac[1][0] = 0.0
    jac[0][1] = 0.0
    jac[1][1] = (2.0/3.0)/(u**(1.0/3.0))
    found_jac = TRUE
  if not found_jac:
    THROW_ARRAY((["unrecognized charts, spacetime=",spacetime,", chart=",chart,", chart2=",chart2]))
  v2[0]=jac[0][0]*v[0]+jac[0][1]*v[1]
  v2[1]=jac[1][0]*v[0]+jac[1][1]*v[1]
  return CLONE_ARRAY_OF_FLOATS(v2) # for python, divorce the new vector from entanglement with components of old

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

