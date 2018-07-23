#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"

xhat=[1.0,0.0,0.0]
yhat=[0.0,1.0,0.0]
zhat=[0.0,0.0,1.0]

def assert_same_vector(u,v):
  test.assert_equal(u[0],v[0])
  

u = euclidean.cross_prod(xhat,yhat)
assert_same_vector(u,zhat)

test.done(verbosity,"euclidean")
