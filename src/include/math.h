#if "LANG" eq "python"
import lambert_w_stuff
from lambert_w_stuff import lambert_w,lambert_w_of_exp
import math_util
from math_util import asinh_of_exp
import math
from math import sin,cos,exp,sinh,cosh,tanh,sqrt,asin,acos,atan2,pi,log,floor
import numpy
import scipy
from numpy import arctanh,arcsinh
from scipy import sign
#endif

#if "LANG" eq "js"
#ifndef IS_BROWSER
#error Symbol IS_BROWSER is not defined in math.h. Make sure to include language.h before math.h.
#endif
#js if (!IS_BROWSER && (typeof Math.karl === 'undefined')) {
#js   /* load() works in rhino, not sure about other engines */
#js   load("lib/math.js");
#js   if (typeof one_over_E === 'undefined') {load("lib/lambertw.js")}
#js   load("lambert_w_stuff.js");
#js   Math.karl = 1;
#js }
#endif


