"""
Low-level math functions for three-dimensional Euclidean geometry.
"""

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "test.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild
#if "LANG" eq "python"
import sys,os,copy
#endif

def rotation_matrix_from_axis_and_angle(theta,u):
  # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  ux = u[0]
  uy = u[1]
  uz = u[2]
  c = cos(theta)
  s = sin(theta)
  return [ \
    [c+ux*ux*(1-c),     ux*uy*(1-c)-uz*s,    ux*uz*(1-c)+uy*s], \
    [uy*ux*(1-c)+uz*s,  c+uy*uy*(1-c),       uy*uz*(1-c)-ux*s], \
    [uz*ux*(1-c)-uy*s,  uz*uy*(1-c)+ux*s,    c+uz*uz*(1-c)] \
  ]

def apply_matrix(m,p):
  return [m[0][0]*p[0]+m[0][1]*p[1]+m[0][2]*p[2],\
          m[1][0]*p[0]+m[1][1]*p[1]+m[1][2]*p[2],\
          m[2][0]*p[0]+m[2][1]*p[1]+m[2][2]*p[2]\
         ]

def dot(a,b):
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def normalize(v):
  n = norm(v)
  v[0] = v[0]/n
  v[1] = v[1]/n
  v[2] = v[2]/n
  return v

def norm(v):
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

def cross_prod(a,b):
  c = [0.0,0.0,0.0]
  c[0] = a[1]*b[2]-a[2]*b[1]
  c[1] = a[2]*b[0]-a[0]*b[2]
  c[2] = a[0]*b[1]-a[1]*b[0]
  return c

