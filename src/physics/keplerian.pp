"""
Low-level routines to compute things for points in the Schwarzschild
spacetime, in my experimental ``Keplerian'' coordinates (t,u,i,j,k),
where u=r^3/2.

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.) 
"""

#include "language.h"
#include "util.h"
#include "math.h"
#include "precision.h"

from io_util import fl

import math_util,angular

def christoffel(p):
  """
  For the Schwarzschild spacetime, compute the Christoffel symbols, in coordinates (a,b,i,j,k).

  The order of indices is as in ctensor:
     symmetric on 1st 2 indices
     contravariant on final index
  See maxima/schwarzschild5.mac. In addition, we put in a fictitious centripetal term that
  keeps us constrained to the unit sphere i^2+j^2+k^2=1.
  """
  return christoffel_raw_maxima_output(p)

def christoffel_raw_maxima_output(p):
  # output of keplerian.mac, plus centripetal terms
  ch = EMPTY3DIM(5)
  u = p[1]
  #------------------------------------------------------
  ch[0][0][1] = 3/(4*u)-3/(4*u**(5/3)) 
  #   ... ^u _t t
  ch[0][1][0] = -u**(1/3)/(3*u**(4/3)-3*u**2) 
  #   ... ^t _t u
  ch[1][0][0] = -u**(1/3)/(3*u**(4/3)-3*u**2) 
  #   ... ^t _u t
  ch[1][1][1] = u**(4/3)/((-9*u)+9*u**(5/3)+u**(1/3)*(3-3*u**2))-(2*u**(2/3))/((-9*u)+9*u**(5/3)+u**(1/3)*(3-3*u**2))+1/((-9*u)+9*u**(5/3)+u**(1/3)*(3-3*u**2)) 
  #   ... ^u _u u
  ch[1][2][2] = 2/(3*u) 
  #   ... ^i _u i
  ch[2][1][2] = 2/(3*u) 
  #   ... ^i _i u
  ch[1][3][3] = 2/(3*u) 
  #   ... ^j _u j
  ch[3][1][3] = 2/(3*u) 
  #   ... ^j _j u
  ch[1][4][4] = 2/(3*u) 
  #   ... ^k _u k
  ch[4][1][4] = 2/(3*u) 
  #   ... ^k _k u
  ch[2][2][1] = (3*u**(1/3))/2-(3*u)/2 
  #   ... ^u _i i
  ch[3][3][1] = (3*u**(1/3))/2-(3*u)/2 
  #   ... ^u _j j
  ch[4][4][1] = (3*u**(1/3))/2-(3*u)/2 
  #   ... ^u _k k
  #------------------------------------------------------
  angular.add_centripetal(ch,p)
  #------------------------------------------------------
  return ch


