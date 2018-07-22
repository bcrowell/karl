#!/usr/bin/python3
#####################

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild
#if "LANG" eq "python"
import sys,os,copy,sqlite3
#endif

verbosity=1

csv=TRUE
csv_file = 'a.csv'

star_catalog = '/usr/share/karl/star_catalog.sqlite'
# Star catalog is built by a script in the directory star_catalog, see README in that directory.

def main():
  r = 9.0
  aberration_table = make_aberration_table(r,1.0e-6,csv,csv_file)
  star_table = make_star_table(star_catalog,aberration_table,r,TRUE,0.0,0.0,2.0)

#--------------------------------------------------------------------------------------------------

def make_aberration_table(r,tol,csv,csv_file):
  """
  Determine a table of optical aberration angles for an observer in the Schwarzschild spacetime.
  Each line of the table is in the format [r,alpha,beta,beta-alpha,f], where r is observer's coordinate
  in units of Schwarzschild radius, alpha is angle as seen by observer if they're in a standard
  state of motion (free fall from rest at infinity), beta is angle on the celestial sphere,
  and f is the factor by which intensities are amplified by lensing. If this is the python version
  running, and csv is true, then we write a copy of the data to csv_file.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  x_obs = [0.0,r,1.0,0.0,0.0]
  # Calculate the velocity vector of the observer, for the standard state of motion:
  v_obs = [1/aa,-sqrt(1-aa),0.0,0.0,0.0]
  #    ... solution of the equations E=Atdot=1, |v|^2=Atdot^2-(1/A)rdot^2=1
  #        We check below that it actually is a solution.
  energy_obs = aa*v_obs[0]
  if abs(energy_obs-1.0)>EPS*10:
    THROW("energy of observer="+str(energy_obs))
  if abs(vector.norm(spacetime,chart,pars,x_obs,v_obs)-1.0)>EPS*10:
    THROW("observer's velocity has bad norm")
  # Calculate the observer's radial vector rho, parallel-transported in from infinity.
  rho = vector.proj(spacetime,chart,pars,x_obs,v_obs,[0.0,1.0,0.0,0.0,0.0])
  # ... take rhat and project out the part orthogonal to the observer's velocity; not yet normalized
  rho = vector.normalize_spacelike(spacetime,chart,pars,x_obs,rho)
  # ... now normalize it to have |rho|^2=-1
  if abs(vector.norm(spacetime,chart,pars,x_obs,rho)+1.0)>EPS*10:
    THROW("norm(rho)!=-1")
  if r>1.0:
    max_le = r/sqrt(aa) # maximum possible L/E for a photon at this r
  else:
    max_le = 10.0
  done = FALSE
  n_angles = 100
  def count_winding(lam,x,v,spacetime,chart,pars):
    if lam==0.0:
      count_winding.winding = 0
      count_winding.last_angle = 0.0
    else:
      angle = atan2(x[3],x[2]) # returns an angle from -pi to pi
      if angle<0.0:
        angle = angle+2.0*MATH_PI
      if count_winding.last_angle>5.0 and angle<1.0:
        count_winding.winding = count_winding.winding+1
      if count_winding.last_angle<1.0 and angle>5.0:
        count_winding.winding = count_winding.winding-1
      count_winding.last_angle = angle
    return count_winding.winding
  count_winding(0.0,[],[],0,0,{})
  table = []
  for in_n_out in range(2): # 0=outward, 1=inward
    for i in range(n_angles):
      z = (float(i)/float(n_angles))
      if in_n_out==1:
        z = 1.0-z
      le = z*max_le # photon's ratio L/E of angular momentum to energy
      # Initial position:
      x = x_obs
      if le==0.0:
        v = [1.0,aa,0.0,0.0,0.0]
      else:
        z = r*r/(le*le)-aa
        if z<0.0:
          break # a photon with angular momentum this big can't exist at this r
        dphi_dr = (1.0/r)/sqrt(z)
        if in_n_out==1:
          dphi_dr = -dphi_dr
        dphi_dt = le*aa/(r*r)
        v = [1.0,dphi_dt/dphi_dr,0.0,dphi_dt,0.0] # tangent vector
      norm = vector.norm(spacetime,chart,pars,x,v)
      if abs(norm)>EPS*10:
        THROW("norm="+str(norm))
      v_perp = vector.proj(spacetime,chart,pars,x_obs,v_obs,v)
      # ... part of photon's velocity orthogonal to observer's velocity
      v_perp = vector.normalize_spacelike(spacetime,chart,pars,x_obs,v_perp)
      # ... normalized
      alpha = acos(-vector.inner_product(spacetime,chart,pars,x_obs,rho,v_perp))
      # ... angle at which the observer says the photon is emitted, see docs
      n = 100
      ndebug=0
      if verbosity>=3:
        ndebug=n/10
      ri = r
      lambda_max = 10.0 # fixme, sort of random
      info = {}
      while True:
        opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,'tol':tol,
              'user_function':count_winding}
        err,final_x,final_v,final_a,final_lambda,info,sigma  = \
                fancy.trajectory_schwarzschild(spacetime,chart,pars,x,v,opt)
        if err!=0:
          if err==RK_INCOMPLETE:
            done = TRUE
            break
          THROW("error: "+str(err))
        rf = final_x[1]
        if IS_NAN(rf):
          THROW("rf is NaN")
        if rf>1.0e9 and rf>100.0*r:
          break
        if rf>10.0 and final_v[1]>0 and rf>ri: # quick and dirty test, far away and getting farther
          f = (rf-ri)/ri # fractional change in r from the iteration we just completed
          if f<0.1:
            lambda_max = lambda_max*(0.1/f) # try to get at least a 10% change in r with each iteration
        if verbosity>=2:
          PRINT("x=",final_x)
        x = final_x
        v = final_v
        ri = rf
      if done:
        break
      w = info['user_data'] # winding number
      beta = atan2(final_x[3],final_x[2]) # returns an angle from -pi to pi
      if beta<0.0:
        beta = beta+2.0*MATH_PI
      beta = beta+w*2.0*MATH_PI
      table.append([r,alpha,beta])
      if verbosity>=1:
        PRINT("r=",r,", alpha=",alpha*180.0/math.pi," deg, beta=",beta*180.0/math.pi," deg, f=",beta-alpha)
    if done:
      break
  table2 = []
  for i in range(len(table)):
    x = CLONE_ARRAY_OF_FLOATS(table[i])
    # x = [r,alpha,beta]
    if i==0:
      ii=0
      jj=1
    else:
      ii=i-1
      jj=i
    alpha = 0.5*(table[ii][1]+table[jj][1]) # average them because we don't want alpha/beta=0/0 at i=0
    beta  = 0.5*(table[ii][2]+table[jj][2])
    dalpha = table[jj][1]-table[ii][1]
    dbeta  = table[jj][2]-table[ii][2]
    f = sin(alpha)*dalpha/(sin(beta)*dbeta) # =dOmega(obs)/dOmega(infinity)=amplification, by Liouville's thm (?)
    x.append(beta-alpha)
    x.append(f)
    table2.append(x)
#if "LANG" eq "python"
  if csv:
    with open(csv_file, 'w') as f:
      for x in table2:
        # x = [r,alpha,beta,beta-alpha,f]
        f.write(",".join(map(lambda u : io_util.fl_n_decimals(u,12), x))+"\n")
    print("table of aberration data written to "+csv_file)
#endif
  return table2

def make_star_table(star_catalog,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag):
  """
  Make a table of stars seen by on observer in the vicinity of a black hole, in a standard state of motion,
  which is free fall from rest at infinity.
  If if_black_hole is false, then the black hole is omitted.
  The parameters ra_out,dec_out are the RA and declination, in radians, of the point on the celestial
  sphere from which the observer has fallen.
  The aberration_table is generated by make_aberration_table(), and contains rows in the format
  [r,alpha,beta,beta-alpha,f]. Here alpha is the observer's angle.
  Max_mag gives the maximum apparent magnitude to include.
  """
  db = sqlite3.connect(star_catalog)
  m = rotation_matrix_observer_to_celestial(ra_out,dec_out,1.0)
  for i in range(len(aberration_table)-1):
    ab1 = aberration_table[i]
    ab2 = aberration_table[i+1]
    alpha1 = ab1[1]
    alpha2 = ab2[1]
    alpha = 0.5*(alpha1+alpha2) # midpoint of strip
    PRINT("alpha=",alpha*180.0/MATH_PI)
    dalpha = alpha2-alpha1
    beta1 = ab1[2]
    beta2 = ab2[2]
    beta = 0.5*(beta1+beta2) # midpoint
    # Break up the 2pi range of azimuthal angles into n_az chunks, so that each chunk within this
    # strip is approximately square.
    n_az = CEIL(2.0*MATH_PI*sin(alpha)/dalpha)
    dphi = 2.0*MATH_PI/n_az
    for j in range(n_az):
      phi1 = dphi*j
      phi2 = dphi*(j+1)
      phi = 0.5*(phi1+phi2)
      # We're now processing a rectangle of the field of view for alpha in [alpha1,alpha2] and phi in [phi1,phi2].
      # Find all four corners of this rectangle, and determine the range of possible RA and dec.
      min_ra = 999.9
      max_ra = 0.0
      min_dec = 999.9
      max_dec = -999.9
      for k in range(2):
        if k==0:
          phi_c = phi1
        else:
          phi_c = phi2
        for l in range(2):
          if l==0:
            beta_c = beta1
          else:
            beta_c = beta2
          ra,dec = beta_to_celestial(beta_c,phi_c,m)
          if ra<min_ra:
            min_ra = ra
          if ra>max_ra:
            max_ra = ra
          if dec<min_dec:
            min_dec = dec
          if dec>max_dec:
            max_dec = dec
      # Do a database search within this range of RA and dec.
      cursor = db.cursor()
      cursor.execute('''SELECT ra,dec,mag,bv FROM stars WHERE mag<=? AND ra>=? AND ra<=? AND dec>=? AND dec<=?''',\
                      (max_mag,min_ra,max_ra,min_dec,max_dec))
      all_rows = cursor.fetchall()
      for row in all_rows:
        ra,dec,mag,bv = row
        PRINT("ra,dec,mag,bv=",ra,dec,mag,bv)
  db.close()

def beta_to_celestial(beta,phi,m):
  """
  An observer is inside the celestial sphere at a certain RA and dec (in radians). They have
  a coordinate system in which beta=0 is the zenith and phi is an azimuthal angle, with phi=0
  being toward the north celestial pole.
  Given a beta and phi, compute the RA and dec of a point below the zenith.
  The matrix m is to have been precomputed using rotation_matrix_observer_to_celestial().
  """
  # Compute the cartesian vector of the point in the observer's frame.
  p = [sin(beta)*cos(phi),sin(beta)*sin(phi),cos(beta)]
  # Rotate to celestial frame:
  p2 =   [m[0][0]*p[0]+m[0][1]*p[1]+m[0][2]*p[2],\
          m[1][0]*p[0]+m[1][1]*p[1]+m[1][2]*p[2],\
          m[2][0]*p[0]+m[2][1]*p[1]+m[2][2]*p[2]\
         ]
  ncp = [0.0,0.0,1.0] # north celestial pole
  cx = [1.0,0.0,0.0] # x axis in celestial coords
  cy = [0.0,1.0,0.0] # x axis in celestial coords
  theta = acos(dot(p2,ncp))
  phi = atan2(dot(p2,cy),dot(p2,cx))
  dec = 0.5*MATH_PI-theta
  ra = phi
  if ra<0.0:
    ra = ra + 2.0*MATH_PI
  return [ra,dec]

def rotation_matrix_observer_to_celestial(ra,dec,direction):
  # See beta_to_celestial() for explanation of what's going on.
  # direction = +-1; if direction is -1, compute the opposite transformation
  # Find the zenith vector, in celestial coords.
  zenith = [cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec)]
  # North celestial pole:
  ncp = [0.0,0.0,1.0]
  # Angle of rotation:
  rot = direction*acos(dot(ncp,zenith))
  # Axis of rotation:
  if rot==0.0:
    axis = ncp # doesn't matter, just avoid division by zero that would otherwise happen in this case
  else:
    axis = normalize(cross_prod(zenith,ncp))
  return rotation_matrix_from_axis_and_angle(rot,axis)

def rotation_matrix_from_axis_and_angle(theta,u):
  # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  ux = u[0]
  uy = u[1]
  uz = u[2]
  c = cos(theta)
  s = sin(theta)
  return [ \
    [c+ux*ux*(1-c),     ux*uy*(1-c)-uz*s,    ux*uz*(1-c)+uy*s], \
    [uy*ux*(1-c)+uz*s,  c+uy*uy*(1-c),       uy*uz*(1-c)-ux*s], \
    [uz*ux*(1-c)-uy*s,  uz*uy*(1-c)+ux*s,    c+uz*uz*(1-c)] \
  ]


def dot(a,b):
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def normalize(v):
  n = norm(v)
  v[0] = v[0]/n
  v[1] = v[1]/n
  v[2] = v[2]/n
  return v

def norm(v):
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

def cross_prod(a,b):
  c = [0.0,0.0,0.0]
  c[0] = a[1]*b[2]-a[2]*b[1]
  c[1] = a[2]*b[0]-a[0]*b[2]
  c[2] = a[0]*b[1]-a[1]*b[0]
  return c

main()
