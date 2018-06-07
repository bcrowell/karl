import io_util
from io_util import strcat

def assert_rel_equal(x,y):
  rel_err = (x-y)/x
  if abs(rel_err)>1.0e-16: raise RuntimeError(strcat(["test failed, x=",x,", y=",y,", rel err=",rel_err]))

