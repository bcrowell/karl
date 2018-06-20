from scipy.special import lambertw

import numpy
import scipy
import math
from numpy import arctanh
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi,log,floor
from scipy import sign

#include "language.h"

###################################################################

#if "LANG" eq "python"
def lambert_w(x):
  """Lambert W function W0."""
  return lambertw(x).real
#endif

###################################################################

def lambert_w_of_exp(u):
  """
  Return W(e^u), where W is the Lambert W function W0.

  The purpose of this function is to handle cases where u is too big
  to allow x=e^u to be stored in floating point. The method used is the
  iterative scheme described in Veberic, https://arxiv.org/abs/1003.1628 , 
  sec. 2.3. The output W of the function satisfies u=ln W+W to within a
  relative error of 10^-16.
  """
  if u<100: return lambert_w(exp(u))
  # Find an initial guess, https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions .
  l2 = log(u) # ln(ln(x))
  nterms = floor(88/l2) # truncate series to avoid underflow; IEEE-754 can represent 2^128; here 88=ln(2^128)
  w = u-l2
  if nterms>=2:
    w += l2/u
    if nterms>=3:
      w += 0.5*l2*(-2+l2)/(u*u)
      if nterms>=4:
        w += l2*(6-9*l2+2*l2*l2)/(6*u*u*u)
  # Iteration:
  # Predetermine how many iterations we need, based on empirical tests.
  if nterms>=5:
    niter=2
  else:
    niter=1
  for i in range(niter):
    z = u-log(w)-w
    q = 2*(1+w)*(1+w+2.0*z/3.0)
    err = (z/(1+w))*((q-z)/(q-2*z))
    w = w*(1+err)
  return w


