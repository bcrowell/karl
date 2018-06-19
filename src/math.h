#if "LANG" eq "python"
import lambert_w
from lambert_w import lambert_w,lambert_w_of_exp
import math_util
from math_util import asinh_of_exp
import math
from math import sin,cos,exp,sinh,cosh,sqrt,asin,acos,atan2,pi,log,floor
import numpy
import scipy
from numpy import arctanh,arcsinh
from scipy import sign
#endif

#if "LANG" eq "js"
if (!IS_BROWSER) {
  /* load() works in rhino, not sure about other engines */
  load("lib/math.js");
  load("lib/lambertw.js");
}
#endif
