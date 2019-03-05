#!/usr/bin/python3

import copy

import numpy as np
import matplotlib.pyplot as plt

import math
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi,log,floor
from numpy import arctanh,arcsinh

from matplotlib import rc

def afuncr(r,b):
  z = r-1.0
  f = z*exp(z)
  a = -arcsinh(math.e*f/sinh(b))
  return a

def afunct(t,b,region):
  # region = 1 for exterior, -1 for interior
  s = sinh(t/2)
  c = cosh(t/2)
  if s==c: return None
  a = region*arcsinh(((s+c)/(s-c))*sinh(b))
  return a

def do_afuncr_lines(scale,n0,maxab,w,which,region,flip_across_main_diag):
  if which=='t':
    n = 3*n0
  else:
    n = n0
  for i in range(n+1):
    u = i*scale # r or t
    b = np.arange(-maxab, maxab, 0.02)
    a = copy.copy(b)
    for j in range(len(b)):
      kill = False
      if which=='r':
        a[j] = afuncr(u,b[j])
      else:
        a[j] = afunct(u,b[j],region)
      kill = kill or abs(a[j])>maxab
      if flip_across_main_diag:
        temp = copy.copy(a[j])
        a[j] = b[j]
        b[j] = temp
      kill = kill or sinh(a[j])*sinh(b[j])>1.0001
      kill = kill or abs(a[j])+abs(b[j])>maxab
      if kill: a[j] = None
    lines = plt.plot(b, a)
    col = 'gray'
    if w>0.75: col='black'
    plt.setp(lines, color=col, linewidth=w)

def do_lines(scale,n,maxab,w,which):
  if which=='r':
    do_afuncr_lines(scale,n,maxab,w,which,0,False)
  else:
    do_afuncr_lines(scale,n,maxab,w,which,-1,False)
    do_afuncr_lines(scale,n,maxab,w,which,1,False)
    do_afuncr_lines(scale,n,maxab,w,which,-1,True)
    do_afuncr_lines(scale,n,maxab,w,which,1,True)

for which in ['r','t']:
  maxab = 5.0   # max value of |a| and |b|
  # --
  n=40
  scale = 0.1
  w = 0.5
  do_lines(scale,n,maxab,w,which)
  # --
  scale = 1.0   # step between values of r
  n = 10        # n+1 values of r
  w = 1.0       # line width in pt?
  do_lines(scale,n,maxab,w,which)

plt.axes().set_aspect('equal')
plt.savefig("graph.svg")
