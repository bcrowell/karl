import io_util
from io_util import strcat

#include "precision.h"

def assert_equal_eps(x,y,eps):
  err = x-y
  if abs(err)>eps: raise RuntimeError(strcat(["test failed, x=",x,", y=",y,", err=",err,", eps=",eps]))

def assert_rel_equal_eps(x,y,eps):
  rel_err = (x-y)/x
  if abs(rel_err)>eps: raise RuntimeError(strcat(["test failed, x=",x,", y=",y,", rel err=",rel_err,
                                                  ", eps=",eps]))

def assert_equal(x,y):
  return assert_equal_eps(x,y,EPS)

def assert_rel_equal(x,y):
  return assert_rel_equal_eps(x,y,EPS)

def done(verbosity,name):
  if verbosity>=1: print("Passed test_"+name)


