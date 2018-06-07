#include "io_util.h"

def assert_equal(x,y):
  err = x-y
  if abs(err)>1.0e-16: raise RuntimeError(strcat(["test failed, x=",x,", y=",y,", err=",err]))

def assert_rel_equal(x,y):
  rel_err = (x-y)/x
  if abs(rel_err)>1.0e-16: raise RuntimeError(strcat(["test failed, x=",x,", y=",y,", rel err=",rel_err]))

def done(verbosity,name):
  if verbosity>=1: print("Passed test_"+name)


