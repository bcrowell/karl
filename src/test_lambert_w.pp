#!/usr/bin/python3

#include "math.h"
#include "init.h"
#include "test.h"

u = 1.0
for i in range(90):
  w = lambert_w_of_exp(u)
  if verbosity>=3: PRINT("====== In test_lambert_w, u=",u,", w=",w)
  assert_rel_equal(u,log(w)+w)
  u = u*2.5 # a little less than e, to try to work out all possible whole-number parts of ln(u)

test.done(verbosity,"lambert_w")
