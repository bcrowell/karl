#!/usr/bin/python3

import lambert_w
from lambert_w import lambert_w,lambert_w_of_exp
import test
from test import assert_rel_equal
import math
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi,log,floor

verbosity=1

u = 1.0
for i in range(90):
  w = lambert_w_of_exp(u)
  if verbosity>=2: print("====== In test_lambert_w, u=",u,", w=",w)
  assert_rel_equal(u,log(w)+w)
  u = u*2.5 # a little less than e, to try to work out all possible whole-number parts of ln(u)

if verbosity>=1: print("Passed test_lambert_w")
