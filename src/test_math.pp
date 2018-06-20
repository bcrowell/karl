#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"

from numpy import arctanh,arcsinh,arccosh

test.assert_rel_equal_eps(1776,cosh(arccosh(1776)),5*EPS)
test.assert_rel_equal_eps(1776,sinh(arcsinh(1776)),5*EPS)
test.assert_rel_equal_eps(0.1776,tanh(arctanh(0.1776)),5*EPS)


done(verbosity,"test_math")
