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
  #------------------------------------------------------
  #------------------------------------------------------
  angular.add_centripetal(ch,p)
  #------------------------------------------------------
  return ch


