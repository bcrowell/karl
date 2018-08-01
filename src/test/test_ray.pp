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

def main():
  r = 30.0/2.0 # the example done in Riazuelo, https://arxiv.org/abs/1511.06025 , p. 7, fig. 1
  test_schwarzschild_standard_observer(r)
  r = 0.5 # test inside horizon
  test_schwarzschild_standard_observer(r) 

def test_schwarzschild_standard_observer(r):
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}

  x_obs,v_obs,rho = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
  aa = 1.0-1.0/r

  # observer's energy is 1:
  assert_equal_eps(aa*v_obs[0],1.0,10*EPS) 

  # rho has norm -1:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,rho),-1.0,10*EPS)

  # v_obs has norm 1:
  assert_equal_eps(vector.norm(spacetime,chart,pars,x_obs,v_obs),1.0,10*EPS)

  # v_obs is orthogonal to rho:
  assert_equal_eps(vector.inner_product(spacetime,chart,pars,x_obs,v_obs,rho),0.0,10*EPS)

  if chart==CH_SCH:
    # rho has a positive r component:
    if r>1.0:
      assert_boolean(rho[1]>0.0,"observer's radial vector has a negative r component in Schwarzschild coordinates")

    # rho has vanishing angular components, observer is supposed to be radially infalling:
    assert_equal_eps(rho[2],0.0,10*EPS)
    assert_equal_eps(rho[3],0.0,10*EPS)
    assert_equal_eps(rho[4],0.0,10*EPS)

    # v_obs is oriented in the future-timelike direction:
    if r>1.0:
      assert_boolean(v_obs[0]>0.0,"observer has r>1 and v_t<=0, i.e., velocity is not future-oriented")
    if r<1.0:
      assert_boolean(v_obs[1]<0.0,"observer has r<1 and v_t>=0, i.e., velocity is not future-oriented")

main()
