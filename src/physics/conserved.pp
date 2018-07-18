import transform,vector

def conserved(chart,x,v):
  """
  For a position and velocity vector given in the Schwarzschild spacetime,
  find the conserved quantities E^-2 (inverse-square energy) and (L/E)^2 (squared
  angular momentum per unit energy).

  The velocity vector does not need to be normalized and may be lightlike or spacelike.
  Returns [e_inv2,le2].
  If r=1 or v_t=0, return [NaN,NaN].

  I'm not sure what is the nicest way to express the L/E so as to retain the three-dimensional
  orientation information and behave correctly and smoothly in all cases, including the case
  of a null velocity vector that becomes lightly spacelike due to rounding errors. The main purpose
  of implementing this is for estimating asymptotics of trajectories that hit the singularity,
  and for that application I don't care about the orientation.
  """
  x2 = transform.transform_point(x,SP_SCH,CH_SCH,{},chart)
  v2 = transform.transform_vector(v,x,SP_SCH,CH_SCH,{},chart)
  r = x2[1]
  aa = 1-1/r
  # Conserved energy is xi^a v_a, where xi^a=(partial_t), so 1/E=(dtau/dt)(1/A), where A=g_tt=1-1/r.
  z = v[0]*aa
  if z==0.0:
    return [NAN,NAN]
  e_inv_sq = vector.norm(SP_SCH,chart,x,v)/(z*z)
  # L/E = r^2 dphi/dt (1/A)
  zz = r*r*pythag(v[2],v[3],v[4])/z # would have exited above if z=0
  return [e_inv_sq,zz*zz]

def pythag(x,y,z):
  return sqrt(x*x+y*y+z*z)
