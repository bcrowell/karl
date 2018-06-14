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
    return metric_sch4(r,sin_theta)
  if (spacetime|chart)==(SP_SCH|CH_AKS): # Schwarzschild metric, in arcsinh-Kruskal coordinates
    return metric_ks4(p)
  raise RuntimeError(strcat(["unrecognized spacetime or chart: ",spacetime," ",chart]))
