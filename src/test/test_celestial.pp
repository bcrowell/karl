#!/usr/bin/python3

#include "language.h"
#include "util.h"
#include "io_util.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "precision.h"

import test,euclidean,celestial

# beta = 90-declination

def test_inverse_rot(ra,dec):
  m_inv = celestial.rotation_matrix_observer_to_celestial(ra,dec,-1.0)
  #print(m_inv)
  beta,phi = celestial.celestial_to_beta(ra,dec,m_inv) # should rotate it to the zenith
  #print("beta=",beta*180.0/MATH_PI," deg, phi=",phi*180.0/MATH_PI," deg")
  test.assert_equal(beta,0.0)

test_inverse_rot(0.0,MATH_PI/2.0)

test_inverse_rot(0.0,0.0)


ra,dec  = celestial.rigel_ra_dec()
test_inverse_rot(ra,dec)

