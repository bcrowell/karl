"""
Low-level routines to compute things for points in the Schwarzschild
spacetime, in Kruskal-Szekeres coordinates, rescaled by an arcsinh,
with the angular space embedded on the unit 2-sphere in a fictitious
additional dimension.  The units are such that the Schwarzschild
radius is 1, i.e., the mass is 1/2. 

Transformations between arcsinh-Kruskal coordinates and Schwarzschild
coordinates, and the associated jacobians, are in transform.pp, except
that transformation from (a,b) to (t,r) are here, in kruskal.aux.

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.) 
"""

#include "util.h"
#include "math.h"
#include "precision.h"

import math_util
from math_util import safe_exp

def metric_ks4(p):
  """
  For the Schwarzschild spacetime, compute the metric.
  Metric is in lower-index form, in +--- signature, with coordinates (a,b,theta,phi).
  """
  g = [[0 for i in range(4)] for j in range(4)]
  t,r,mu = aux(p[0],p[1])
  g[0][1] = mu
  g[1][0] = mu
  r2 = r*r
  sin_theta = sin(p[2])
  g[2][2] = -r2
  g[3][3] = -r2*sin_theta*sin_theta
  return g

def aux(a,b):
  """
  Compute Schwarzschild t and r, and metric element mu, for a point given in rescaled Kruskal coordinates (a,b).
  Return [t,r,mu].

  Although r is an analytic function of the coordinates throughout the entire
  maximal extension of the spacetime, numerical issues force us to do some case-splitting.
  The quantity mu is the coefficient in the line element, ds^2 = 2 mu da db -...,
  i.e., it's the matrix element g_ab of the metric. Returns t=None for points on a horizon.
  Returns None for all three variables if a and b lie beyond the singularity.
  """
  if a<=0: # Flip to positive a in order to simplify some later computations.
    return aux(-a,-b)
  # From now on, we know a>=0.
  e2a = safe_exp(-2*a)
  e2b = safe_exp(-2*abs(b))
  if b>=0.0: # Can I simplify the logic by eliminating this?
    # interior (region II) or horizon of II
    r,mu = small_r_mu(a,b)
  else:
    # exterior, region I
    # From now on, we know a>0 and b<0.
    # Here, if r is even moderately large the coordinates v and w can be too big to store in
    # floating point, so we use a different algorithm.
    f = (1.0-e2a)*(1.0-e2b)
    u = a-b+log(f/4.0)-1.0
    r = 1.0+lambert_w_of_exp(u)
    # Compute mu:
    mu = (1.0+e2a)*(1.0+e2b)*(1/(2*math.e*r))*exp(a-b-(r-1))
  # Compute t:
  if a!=0 and b!=0 and r is not None:
    t = a-abs(b)+log((1-e2a)/(1-e2b))
  else:
    t = None
  return [t,r,mu]

def small_r_mu(a,b):
  """
  Returns [r,mu]. Works only at small r and t, produces overflows elsewhere. Returns [None,None] in
  regions beyond the singularity.
  """
  v=sinh(a); w=sinh(b)
  rho = -v*w
  if rho<-1.0: return [None,None]
  r = 1+lambert_w(rho/math.e)
  mu = 2.0*((r-1)/(r*rho))*cosh(a)*cosh(b)
  return [r,mu]

def metric(p):
  """
  For the Schwarzschild spacetime, compute the metric.
  Metric is in lower-index form, in +--- signature, with coordinates (a,b,i,j,k).
  """
  g = [[0 for i in range(5)] for j in range(5)]
  t,r,mu = aux(p[0],p[1])
  g[0][1] = mu
  g[1][0] = mu
  r2 = r*r
  g[2][2] = -r2
  g[3][3] = -r2
  g[4][4] = -r2
  return g

def christoffel(p):
  """
  For the Schwarzschild spacetime, compute the Christoffel symbols, in coordinates (a,b,i,j,k).

  The order of indices is as in ctensor:
     symmetric on 1st 2 indices
     contravariant on final index
  See maxima/schwarzschild5.mac. In addition, we put in a fictitious centripetal term that
  keeps us constrained to the unit sphere i^2+j^2+k^2=1.
  """
  #return christoffel_raw_maxima_output(p)
  return christoffel_massaged_maxima_output(p)

def christoffel_massaged_maxima_output(p):
  a=p[0]; b=p[1]
  # Based on christoffel_raw_maxima_output(), but simplified and massaged into a form that won't cause overflows.
  # FIXME: haven't fixed overflows yet.
  if a<0: # regions III and IV
    p2 = copy.copy(p)
    p2[0] = -a
    p2[1] = -b
    return flip_christoffel(christoffel_massaged_maxima_output(p2))
  # From now on, we know we're in region I or II.
  ch = [[[0 for i in range(5)] for j in range(5)] for k in range(5)]
  #------------------------------------------------------
  t,r,mu = aux(a,b)
  r2 = r*r
  r2e = (1/r2)*exp(-r)
  #------------------------------------------------------
  ch[0][0][0] = math_util.safe_tanh(a)+(1.0/r+1.0/r2)*exp(-r)*cosh(a)*sinh(b)
  ch[1][1][1] = math_util.safe_tanh(b)+(1.0/r+1.0/r2)*exp(-r)*cosh(b)*sinh(a)
  #--
  z = -r2e*cosh(a)*sinh(b)
  ch[0][2][2] = z
  ch[2][0][2] = z
  ch[0][3][3] = z
  ch[3][0][3] = z
  ch[0][4][4] = z
  ch[4][0][4] = z
  #--
  z = -r2e*cosh(b)*sinh(a)
  ch[1][2][2] = z
  ch[2][1][2] = z
  ch[1][3][3] = z
  ch[3][1][3] = z
  ch[1][4][4] = z
  ch[4][1][4] = z
  #--
  z = -0.5*r*math_util.safe_tanh(a)
  ch[2][2][0] = z
  ch[3][3][0] = z
  ch[4][4][0] = z
  #--
  z = -0.5*r*math_util.safe_tanh(b)
  ch[2][2][1] = z
  ch[3][3][1] = z
  ch[4][4][1] = z
  #------------------------------------------------------
  add_centripetal(ch,p)
  #------------------------------------------------------
  return ch  

def flip_christoffel(ch):
  """
  Given the Christoffel symbols at (a,b), change them, in place, to the ones at (-a,-b), and return the result.
  """
  # For every index that equals a or b, we get one sign flip.
  for i in range(5):
    for j in range(5):
      for k in range(5):
        ch[i][j][k] = flip_a_christoffel_helper(i,j,k)*ch[i][j][k]
  return ch

def flip_a_christoffel_helper(i,j,k):
  s = 1.0
  if i<=1: s = -s
  if j<=1: s = -s
  if k<=1: s = -s
  return s

def christoffel_raw_maxima_output(p):
  a=p[0]; b=p[1]
  # output of kruskal5.mac, plus centripetal terms
  ch = [[[0 for i in range(5)] for j in range(5)] for k in range(5)]
  #------------------------------------------------------
  ch[0][0][0] = (sinh(a)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2+(cosh(a)**2*sinh(b)+2*sinh(a)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1))*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(a)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*cosh(a)**2*sinh(b))/(cosh(a)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*cosh(a)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+cosh(a)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^a _a a
  ch[0][2][2] = -(cosh(a)*sinh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^i _a i
  ch[2][0][2] = -(cosh(a)*sinh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^i _i a
  ch[0][3][3] = -(cosh(a)*sinh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^j _a j
  ch[3][0][3] = -(cosh(a)*sinh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^j _j a
  ch[0][4][4] = -(cosh(a)*sinh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^k _a k
  ch[4][0][4] = -(cosh(a)*sinh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^k _k a
  ch[1][1][1] = (sinh(b)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2+(sinh(a)*cosh(b)**2+2*sinh(b)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1))*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(b)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*sinh(a)*cosh(b)**2)/(cosh(b)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*cosh(b)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+cosh(b)*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^b _b b
  ch[1][2][2] = -(sinh(a)*cosh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^i _b i
  ch[2][1][2] = -(sinh(a)*cosh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^i _i b
  ch[1][3][3] = -(sinh(a)*cosh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^j _b j
  ch[3][1][3] = -(sinh(a)*cosh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^j _j b
  ch[1][4][4] = -(sinh(a)*cosh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^k _b k
  ch[4][1][4] = -(sinh(a)*cosh(b))/(exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)+2*exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+exp(lambert_w(-exp(-1)*sinh(a)*sinh(b))+1)*lambert_w(-exp(-1)*sinh(a)*sinh(b))**2) 
  #   ... ^k _k b
  ch[2][2][0] = -(sinh(a)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(a))/(2*cosh(a)) 
  #   ... ^a _i i
  ch[2][2][1] = -(sinh(b)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(b))/(2*cosh(b)) 
  #   ... ^b _i i
  ch[3][3][0] = -(sinh(a)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(a))/(2*cosh(a)) 
  #   ... ^a _j j
  ch[3][3][1] = -(sinh(b)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(b))/(2*cosh(b)) 
  #   ... ^b _j j
  ch[4][4][0] = -(sinh(a)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(a))/(2*cosh(a)) 
  #   ... ^a _k k
  ch[4][4][1] = -(sinh(b)*lambert_w(-exp(-1)*sinh(a)*sinh(b))+sinh(b))/(2*cosh(b)) 
  #   ... ^b _k k
  #------------------------------------------------------
  add_centripetal(ch,p)
  #------------------------------------------------------
  return ch

def christoffel_raw_maxima_output_first_try(p):
  a=p[0]; b=p[1]
  # output of kruskal5.mac, plus centripetal terms
  ch = [[[0 for i in range(5)] for j in range(5)] for k in range(5)]
  #------------------------------------------------------
  ch[0][0][0] = (sinh(a)*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2+(2*sinh(a)*exp(lambert_w(exp((-b)+a-1))+b+1)-exp(a)*cosh(a))*lambert_w(exp((-b)+a-1))+sinh(a)*exp(lambert_w(exp((-b)+a-1))+b+1)-2*exp(a)*cosh(a))/(cosh(a)*exp(lambert_w(exp((-b)+a-1))+b+1)+2*cosh(a)*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+cosh(a)*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^a _a a
  ch[0][2][2] = exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^i _a i
  ch[2][0][2] = exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^i _i a
  ch[0][3][3] = exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^j _a j
  ch[3][0][3] = exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^j _j a
  ch[0][4][4] = exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^k _a k
  ch[4][0][4] = exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^k _k a
  ch[1][1][1] = (exp(lambert_w(exp((-b)+a-1))+b+1)*sinh(b)*lambert_w(exp((-b)+a-1))**2+(exp(a)*cosh(b)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*sinh(b))*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*sinh(b)+2*exp(a)*cosh(b))/(exp(lambert_w(exp((-b)+a-1))+b+1)*cosh(b)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*cosh(b)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*cosh(b)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^b _b b
  ch[1][2][2] = -exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^i _b i
  ch[2][1][2] = -exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^i _i b
  ch[1][3][3] = -exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^j _b j
  ch[3][1][3] = -exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^j _j b
  ch[1][4][4] = -exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^k _b k
  ch[4][1][4] = -exp(a)/(exp(lambert_w(exp((-b)+a-1))+b+1)+2*exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))+exp(lambert_w(exp((-b)+a-1))+b+1)*lambert_w(exp((-b)+a-1))**2) 
  #   ... ^k _k b
  ch[2][2][0] = -(exp(-b)*(exp(a)*lambert_w(exp((-b)+a-1))+exp(a)))/(2*cosh(a)*cosh(b)) 
  #   ... ^a _i i
  ch[2][2][1] = (exp(-b)*(exp(a)*lambert_w(exp((-b)+a-1))+exp(a)))/(2*cosh(a)*cosh(b)) 
  #   ... ^b _i i
  ch[3][3][0] = -(exp(-b)*(exp(a)*lambert_w(exp((-b)+a-1))+exp(a)))/(2*cosh(a)*cosh(b)) 
  #   ... ^a _j j
  ch[3][3][1] = (exp(-b)*(exp(a)*lambert_w(exp((-b)+a-1))+exp(a)))/(2*cosh(a)*cosh(b)) 
  #   ... ^b _j j
  ch[4][4][0] = -(exp(-b)*(exp(a)*lambert_w(exp((-b)+a-1))+exp(a)))/(2*cosh(a)*cosh(b)) 
  #   ... ^a _k k
  ch[4][4][1] = (exp(-b)*(exp(a)*lambert_w(exp((-b)+a-1))+exp(a)))/(2*cosh(a)*cosh(b)) 
  #   ... ^b _k k
  #------------------------------------------------------
  add_centripetal(ch,p)
  #------------------------------------------------------
  return ch

def add_centripetal(ch,p):
  """
  Modifies ch in place by adding the centripetal parts.
  """
  i=p[2] ; j=p[3]; k=p[4]
  xi2 = i*i+j*j+k*k; # should normally be very close to 1
  for m in range(2,5): # upper index
    z = p[m]
    for n in range(2,5): # lower indices
      ch[n][n][m] += z/xi2

def christoffel4(p):
  """
  For the Schwarzschild spacetime, compute the Christoffel symbols, in coordinates (a,b,theta,phi).

  The order of indices is as in ctensor:
     symmetric on 1st 2 indices
     contravariant on final index
  """
  ch = [[[0 for i in range(4)] for j in range(4)] for k in range(4)]
  a=p[0] ; b=p[1]
  t,r,mu = aux(a,b)
  ta = math_util.safe_tanh(a)
  tb = math_util.safe_tanh(b)
  # FIXME: In the following, we'll get overflows if a or b is large; massage the equations to avoid this problem.
  ca=cosh(a); cb=cosh(b); sa=sinh(a); sb=sinh(b)
  q = (1.0/r+1.0/(r*r))
  big_b = (4.0/r)*exp(-r)
  eighth_r = r/8.0
  sin_theta = sin(theta)
  sin2_theta = sin_theta**2
  cos_theta = cos(theta)
  ch[0][0][0] = ta+q*exp(-r)*ca*sb      # ^a_aa
  ch[1][1][1] = tb+q*exp(-r)*cb*sa      # ^b_bb
  ch[0][2][2] = -(big_b/(4*r))*sb/ca    # ^theta_a theta
  ch[2][0][2] = ch[0][2][2]             # ^theta_theta a
  ch[0][3][3] = ch[0][2][2]             # ^phi_a phi
  ch[3][0][3] = ch[0][2][2]             # ^phi_phi a
  ch[1][2][2] = -(big_b/(4*r))*sa/cb    # ^theta_b theta
  ch[1][3][3] = ch[1][2][2]             # ^phi_b phi
  ch[2][2][0] = -eighth_r*ta            # ^a_theta theta
  ch[2][2][1] = -eighth_r*tb            # ^a_theta theta
  ch[3][3][0] = -eighth_r*ta*sin2_theta   # ^a_theta theta
  ch[3][3][1] = -eighth_r*tb*sin2_theta   # ^a_theta theta
  ch[3][3][2] = -sin_theta*cos_theta    # ^theta_phi phi
  ch[2][3][3] = cos_theta/sin_theta     # ^phi_theta phi
  ch[3][2][3] = ch[2][3][3]             # ^phi_phi theta
  return ch


