import math
from math import sin,cos,tan,exp,sinh,cosh,sqrt,asin,acos,tanh,atan2,pi,log,floor
import numpy
import scipy
from numpy import arctanh,arcsinh,arccosh
from scipy import sign

#include "language.h"
#include "precision.h"

def safe_exp(x):
  """
  Return exp(x), or 0 if x is so large and negative that exp(x) is negligible compared to unity,
  and might have caused an underflow.
  """
  if x<LN_EPS:
    return 0.0
  return exp(x)

def safe_tanh(x):
  """
  Return tanh(x), without overflows, underflows, or unnecessarily inefficient evaluation of negligible stuff.
  """
  if x<0.0:
    return -safe_tanh(-x)
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
  if u>LN_EPS:
    return y # correction term will be negligible compared to unity
  g = exp(-2.0*u) # 1/x^2, where x is the input to the arcsinh
  m = 1.0
  coeffs = [1.0/4.0,-3.0/32.0,15.0/288.0,-105.0/3072.0]
  for n in range(1,6):
    m = m*g
    term = coeffs[n-1]*m
    # print("calculating correction, n=",n,", term=",term)
    if term<EPS:
      return y
    y = y + term
  return y # If we get here, the approx. may not be very good.

def force_into_range(x,a,b):
  if x<a:
    return a
  if x>b:
    return b
  return x

def linear_interp(x1,x2,y1,y2,x):
  # Consider using the _safe version of this, below.
  return ((y2-y1)/(x2-x1))*(x-x1)+y1

def linear_interp_safe(x1,x2,y1,y2,x,allow_extrap,sanity_check_extrap,max_extrap):
  unsafe = FALSE
  requires_extrap = not (x1<=x and x2>=x)
  if requires_extrap:
    if allow_extrap:
      if sanity_check_extrap:
        if x<x1:
          unsafe = abs(x-x1)>max_extrap
        else:
          unsafe = abs(x-x2)>max_extrap
    else:
      unsafe = TRUE
  result = linear_interp(x1,x2,y1,y2,x)
  return [result,unsafe]

def linear_interp_from_table_safe(table,x_col,y_col,x,allow_extrap,sanity_check_extrap,max_extrap):
  """
  Returns [result,unsafe].
  """
  return linear_interp_from_table_recurse(table,x_col,y_col,x,0,len(table)-1,\
              allow_extrap,sanity_check_extrap,max_extrap)

def linear_interp_from_table(table,x_col,y_col,x):
  # Better to use the _safe version of this, above.
  result,unsafe = linear_interp_from_table_recurse(table,x_col,y_col,x,0,len(table)-1,TRUE,FALSE,0.0)
  return result

def linear_interp_from_table_recurse(table,x_col,y_col,x,i,j,allow_extrap,sanity_check_extrap,max_extrap):
  # We don't normally call this directly, we call linear_interp_from_table() or linear_interp_from_table_safe(),
  # which calls this routine.
  # Do a binary search through the table, so this is O(log(n)).
  # During recursion, the i,j parameters keep track of what range of row numbers we've narrowed in on.
  # The x column has to be sorted in ascending order.
  # Returns [result,unsafe], where unsafe is a boolean that tells us whether extrapolation was necessary
  # and violated the sanity limits.
  if j==i+1:
    unsafe = FALSE
    requires_extrap = not (table[i][x_col]<=x and table[j][x_col]>=x)
    if requires_extrap:
      if allow_extrap:
        if sanity_check_extrap:
          if x<table[i][x_col]:
            unsafe = abs(x-table[i][x_col])>max_extrap
          else:
            unsafe = abs(x-table[j][x_col])>max_extrap
      else:
        unsafe = TRUE
    result = linear_interp(table[i][x_col],table[j][x_col],table[i][y_col],table[j][y_col],x)
    return [result,unsafe]
  k=(i+j)//2
  if table[k][x_col]>x:
    return linear_interp_from_table_recurse(table,x_col,y_col,x,i,k,allow_extrap,sanity_check_extrap,max_extrap)
  else:
    return linear_interp_from_table_recurse(table,x_col,y_col,x,k,j,allow_extrap,sanity_check_extrap,max_extrap)

