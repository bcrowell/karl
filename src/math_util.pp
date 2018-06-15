import math
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi,log,floor
import numpy
import scipy
from numpy import arctanh,arcsinh
from scipy import sign

#include "precision.h"

def safe_exp(x):
  """
  Return exp(x), or 0 if x is so large and negative that exp(x) is negligible compared to unity,
  and might have caused an underflow.
  """
  if x<LN_EPS: return 0.0
  return exp(x)

def safe_tanh(x):
  """
  Return tanh(x), without overflows, underflows, or unnecessarily inefficient evaluation of negligible stuff.
  """
  if x<0.0: return -safe_tanh(-x)
  # From here on, we know x>=0.
  z = safe_exp(-2*x)
  return (1.0-z)/(1.0+z)

def asinh_of_exp(u):
  """
  Compute asinh(e^u), using asymptotics if u is large.

  Series expansion: https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Series_expansions
  """
  if u<0.5*LN_MAX_FLOAT:
    # ... IEEE-754 floating point can store exp(88); asinh probably uses the square of the input internally
    return arcsinh(exp(u))
  y = log(2.0)+u
  if u>LN_EPS: return y # correction term will be negligible compared to unity
  g = exp(-2.0*u) # 1/x^2, where x is the input to the arcsinh
  m = 1.0
  coeffs = [1.0/4.0,-3.0/32.0,15.0/288.0,-105.0/3072.0]
  for n in range(1,6):
    m = m*g
    term = coeffs[n-1]*m
    # print("calculating correction, n=",n,", term=",term)
    if term<EPS: return y
    y = y + term
  return y # If we get here, the approx. may not be very good.

