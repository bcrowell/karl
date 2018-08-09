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

  The relevant equations are written up more legibly in the docs.
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

def sch_is_in_future_light_cone(x,v):
  """
  Given vector v in Schwarzschild coordinates, determines whether it's in the future light cone.
  Returns [is_future,a_plus_b], where is_future is a boolean, and ab_sum is the value of v_a+v_b in
  KS coordinates. The boolean can be unreliable due to rounding if v is lightline; if you have a vector
  that is known lightlike and you just want to know whether it's future-oriented, it's better to test
  whether the sign of ab_sum is positive.
  """
  # There is probably a more efficient way to do this, but this way seems bulletproof and manifestly correct.
  v_kruskal = transform_vector(v,x,SP_SCH,CH_SCH,{},CH_AKS)
  va = v_kruskal[0]
  vb = v_kruskal[1]
  tol = 10*EPS
  return [(va>=-tol and vb>=-tol),va+vb]

def kruskal_to_time_zero(x,v,spacetime_or_chart,force):
  """
  Change kruskal coordinates to equivalent coordinates that correspond to Schwarzschild time t=0.
  If boolean input force is false, then nothing is done unless the point is the type of point near
  the horizon or photon sphere at large t for which this operation is likely to be helpful to precision.
  If the point is on the horizon, this operation is a no-op. As a convenience, this function can
  also be called when coords are not AKS, and then it's also a no-op.
  Returns [x2,v2,dt,did_it], where dt is the change in the Schwarzschild time that resulted from the
  transformation, and dit_it is true if anything actually happened.
  """
  if spacetime_or_chart!=(SP_SCH|CH_AKS):
    return [x,v,0.0,FALSE]
  a = x[0]
  b = x[1]
  # --
  big = 100.0
  # ...If this parameter is very big, like 10^6, then this code never gets executed and can't serve its
  #    purpose, while if I made it very small, like 1, then this code would get executed frequently,
  #    causing more rounding errors.
  #--
  do_it = TRUE
  if not force:
    do_it = (abs(a)>big*abs(b)) or (abs(b)>big*abs(a))
  if not do_it:
    return [x,v,0.0,FALSE]
  if b==0.0 or a==0.0: # can't define this operation for points on the horizon
    return [x,v,0.0,FALSE]
  # If we're in region III or IV, then simply doing AKS->Schwarzschild->AKS would lose that information.
  # In that case, transform into I or II and recurse.
  if a<0.0:
    # Invert in the origin in AKS.
    x2 = CLONE_ARRAY_OF_FLOATS(x)
    v2 = CLONE_ARRAY_OF_FLOATS(v)
    x2[0] = -x2[0]
    x2[1] = -x2[1]
    v2[0] = -v2[0]
    v2[1] = -v2[1]
    x3,v3,dt,did_it = kruskal_to_time_zero(x2,v2,spacetime_or_chart,force)
    x3[0] = -x3[0]
    x3[1] = -x3[1]
    v3[0] = -v3[0]
    v3[1] = -v3[1]
    return [x3,v3,dt,did_it]
  # Beyond this point, we're guaranteed to be in region I or II, and not at a horizon.
  # The following is a composition of three transformations: (1) AKS to Schwarzschild, (2) a time
  # translation, and (3) back to AKS. The jacobian for 2 is just the identity matrix, so the velocity
  # is not affected by it.
  x_s = transform_point(x,SP_SCH,CH_AKS,{},CH_SCH)
  v_s = transform_vector(v,x,SP_SCH,CH_AKS,{},CH_SCH)
  dt = x_s[0]
  x_s[0] = 0.0
  x2 = transform_point(x_s,SP_SCH,CH_SCH,{},CH_AKS)
  v2 = transform_vector(v_s,x_s,SP_SCH,CH_SCH,{},CH_AKS)
  return [x2,v2,dt,TRUE]
