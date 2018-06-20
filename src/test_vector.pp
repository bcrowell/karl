#!/usr/bin/python3

#include "language.h"
#include "math.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "precision.h"

import vector,transform

def test_norm_schwarzschild_vs_kruskal(t,r,dt,dr,eps):
  p_s = [t,r,0,0] # point in Schwarzschild coords
  norm_s = vector.norm4(SP_SCH,CH_SCH,p_s,[dt,dr,0,0])
  # Transform to arcsinh-Kruskal and compare norms:
  a,b = transform.schwarzschild_to_kruskal(t,r)
  a2,b2 = transform.schwarzschild_to_kruskal(t+dt,r+dr)
  p_k = [a,b,0,0]
  norm_k = vector.norm4(SP_SCH,CH_AKS,p_k,[a2-a,b2-b,0,0])
  if verbosity>=3: print("test_norm_schwarzschild_vs_kruskal: t=",t,", r=",r,", a=",a,", b=",b,
            ", norm_s=",norm_s,", norm_k=",norm_k)
  test.assert_rel_equal_eps(norm_s,norm_k,eps)

d = sqrt(EPS)
max_err = 20*d
test_norm_schwarzschild_vs_kruskal(0.111,2.0,d,0.777*d,max_err) # random point and displacement in region I
test_norm_schwarzschild_vs_kruskal(0.111,0.5,d,0.777*d,max_err) # ... region II

test.done(verbosity,"vector")

