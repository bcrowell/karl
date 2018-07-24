#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"

import euclidean

xhat=[1.0,0.0,0.0]
yhat=[0.0,1.0,0.0]
zhat=[0.0,0.0,1.0]

ident = [xhat,yhat,zhat] # identity matrix

random_unit = euclidean.normalize([1.0,3.4567,1.776])

rot90x = euclidean.rotation_matrix_from_axis_and_angle(0.5*MATH_PI,xhat)
rot180x = euclidean.rotation_matrix_from_axis_and_angle(MATH_PI,xhat)
rot1radiany = euclidean.rotation_matrix_from_axis_and_angle(1.0,yhat)

def assert_same_vector(u,v):
  test.assert_equal(u[0],v[0])
  test.assert_equal(u[1],v[1])
  test.assert_equal(u[2],v[2])
  
assert_same_vector(euclidean.cross_prod(xhat,yhat),zhat)
assert_same_vector(euclidean.cross_prod(yhat,zhat),xhat)
assert_same_vector(euclidean.cross_prod(zhat,xhat),yhat)

test.assert_equal(euclidean.determinant(ident),1.0)
test.assert_equal(euclidean.determinant(rot90x),1.0)
test.assert_equal(euclidean.determinant(euclidean.rotation_matrix_from_axis_and_angle(0.123,random_unit)),1.0)

assert_same_vector(euclidean.apply_matrix(ident,xhat),xhat)
assert_same_vector(euclidean.apply_matrix(rot90x,yhat),zhat)
test.assert_equal(euclidean.dot(xhat,euclidean.apply_matrix(rot1radiany,xhat)),cos(1.0))
test.assert_equal(euclidean.dot(yhat,euclidean.apply_matrix(rot180x,yhat)),-1.0)

test.done(verbosity,"euclidean")
