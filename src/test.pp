import io_util
from io_util import strcat
#js load("io_util.js");

#include "language.h"
#include "math.h"
#include "precision.h"

def assert_equal_eps(x,y,eps):
  err = x-y
  if IS_NAN(x) or IS_NAN(y) or abs(err)>eps:
    THROW(io_util.strcat(["test failed, x=",x,", y=",y,", err=",err,", eps=",eps]))

def assert_rel_equal_eps(x,y,eps):
  if x==0.0 and y==0.0:
    return
  if x==0.0:
    return assert_rel_equal_eps(y,x,eps) # avoid division by zero
  rel_err = (x-y)/x
  if IS_NAN(x) or IS_NAN(y) or abs(rel_err)>eps:
    THROW(io_util.strcat(["test failed, x=",x,", y=",y,", rel err=",rel_err,", eps=",eps]))

def assert_equal(x,y):
  return assert_equal_eps(x,y,2.0*EPS)

def assert_rel_equal(x,y):
  return assert_rel_equal_eps(x,y,2.0*EPS)

def assert_rel_equal_eps_vector(x,y,eps):
  for i in range(len(x)):
    assert_rel_equal_eps(x[i],y[i],eps)

def done(verbosity,name):
  if verbosity>=1:
    print("Passed test_"+name)

verbosity=1
