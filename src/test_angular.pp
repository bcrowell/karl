#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "init.h"
#include "test.h"

import angular

x = [0,0,1,2,3]
n = angular.renormalize(x)
test.assert_equal(n[2]**2+n[3]**2+n[4]**2,1)

x = [0,0,1,0,0]
v = [0,0,3,700,0]
v = angular.make_tangent(x,v)
test.assert_equal(v[2],0)
test.assert_equal(v[3],700)

test.done(verbosity,"angular")
