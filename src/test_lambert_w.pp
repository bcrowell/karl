#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "precision.h"
#include "init.h"
#include "test.h"

import lambert_w_stuff

#if "LANG" eq "python"
import ctypes
import c_libs
#...don't import karl_c_lib this way, because you get a copy of it in whatever state it was in?
#    https://stackoverflow.com/a/142601/1142217 -- comment by Paul Whipp
#endif

z = 0.1492
test.assert_rel_equal_eps(z,lambert_w(z*exp(z)),5*EPS)
test.assert_rel_equal_eps(z,c_libs.karl_c_lib.veberic_lambert_w(z*exp(z)),5*EPS)

u = 1.0
for i in range(90):
  w = lambert_w_stuff.lambert_w_of_exp(u)
  if verbosity>=3:
    PRINT("====== In test_lambert_w, u=",u,", w=",w)
  assert_rel_equal(u,log(w)+w)
  u = u*2.5 # a little less than e, to try to work out all possible whole-number parts of ln(u)

test.done(verbosity,"lambert_w")
