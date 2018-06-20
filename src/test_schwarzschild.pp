#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "init.h"
#include "test.h"

import schwarzschild

def test_christoffel_symmetry(ch):
  l = len(ch)
  for i in range(l):
    for j in range(l):
      for k in range(l):
        test.assert_equal(ch[i][j][k],ch[j][i][k])

# Exercise as many lines of code as possible, using random coordinates.
t = 3.3
r = 2.2
sin_theta = 1/sqrt(2.0)
cos_theta = 1/sqrt(2.0)
# A set of i,j,k that lie on the unit sphere:
i=0.6
j=0.8
k=0.0

schwarzschild.metric_sch4(r,sin_theta)
schwarzschild.sigma(r)
ch4 = schwarzschild.christoffel_sch4(t,r,sin_theta,cos_theta)
ch5 = schwarzschild.christoffel([t,r,i,j,k])

# test symmetry of Christoffel symbols
test_christoffel_symmetry(ch4)
test_christoffel_symmetry(ch5)

# As far as testing whether the results are actually correct, see
# tests in test_runge_kutta.

test.done(verbosity,"schwarzschild")

