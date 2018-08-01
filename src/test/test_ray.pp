#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"
#include "spacetimes.h"

import test,vector,ray,star_properties

spacetime = SP_SCH
chart = CH_SCH
pars = {}
r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1

x_obs,v_obs,rho = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
aa = 1.0-1.0/r
assert_equal_eps(aa*v_obs[0],1.0,10*EPS) # observer's energy is 1
assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,rho),-1.0,10*EPS) # rho has norm -1
assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,v_obs),1.0,10*EPS) # v_obs has norm 1

