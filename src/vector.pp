#include "language.h"
#include "util.h"
#include "math.h"
#include "spacetimes.h"

import schwarzschild,angular,kruskal

def norm4(spacetime,chart,p,v):
  """
  Returns the norm of the vector, in the standard 4-dimensional coordinates.
  """
  g = get_metric4(spacetime,chart,p)
  n = 0.0
  for i in range(4):
    for j in range(4):
      n += g[i][j]*v[i]*v[j]
  return n

def get_metric4(spacetime,chart,p):
  """
  Returns the lower-index form of the metric at the point p, in the standard 4-dimensional
  coordinates.
  """
  if (spacetime|chart)==(SP_SCH|CH_SCH): # Schwarzschild metric, in Schwarzschild coordinates
    r = p[1]
    sin_theta = sin(p[2])
    return schwarzschild.metric_sch4(r,sin_theta)
  if (spacetime|chart)==(SP_SCH|CH_AKS): # Schwarzschild metric, in arcsinh-Kruskal coordinates
    return kruskal.metric_ks4(p)
  raise RuntimeError(io_util.strcat(["unrecognized spacetime or chart: ",spacetime," ",chart]))

def scalar_mult(v0,s):
  """
  Returns a copy of v, which has been multiplied by the scalar s. Works for 4 or 5 dimensions.
  """
  v = copy.copy(v0)
  for i in range(len(v)):
    v[i] = v[i]*s
  return v

def normalize(spacetime,chart,p,v):
  """
  Returns a copy of v, which has been normalized. Works for 4 or 5 dimensions.
  The vector has to be timelike.
  """
  n = norm(spacetime,chart,p,v)
  s = 1/sqrt(n)
  return scalar_mult(v,s)

def norm(spacetime,chart,p,v):
  """
  Returns the norm of the vector, in 5-dimensional coordinates.
  """
  g = get_metric(spacetime,chart,p)
  n = 0.0
  for i in range(5):
    for j in range(5):
      n += g[i][j]*v[i]*v[j]
  return n

def get_metric(spacetime,chart,p):
  """
  Returns the lower-index form of the metric at the point p, in 5-dimensional
  coordinates.
  """
  if (spacetime|chart)==(SP_SCH|CH_SCH): # Schwarzschild metric, in Schwarzschild coordinates
    r = p[1]
    return schwarzschild.metric(r)
  if (spacetime|chart)==(SP_SCH|CH_AKS): # Schwarzschild metric, in arcsinh-Kruskal coordinates
    return kruskal.metric(p)
  raise RuntimeError(io_util.strcat(["unrecognized spacetime or chart: ",spacetime," ",chart]))
