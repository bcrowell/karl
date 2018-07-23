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
import sys,os,copy,sqlite3,csv
#endif

verbosity=1

write_csv=TRUE
csv_file = 'a.csv'

#star_catalog = '/usr/share/karl/star_catalog.sqlite'
star_catalog = '/usr/share/karl/mag7.sqlite'
# Star catalog is built by a script in the directory star_catalog, see README in that directory.

def main():
  r = 100.0
  # falling inward from the direction of Rigel, https://en.wikipedia.org/wiki/Rigel :
  ra_out,dec_out = rigel_ra_dec()
  aberration_table = make_aberration_table(r,1.0e-6,write_csv,csv_file)
  star_table = make_star_table(star_catalog,aberration_table,r,TRUE,ra_out,dec_out,5.0,write_csv,csv_file)

#--------------------------------------------------------------------------------------------------

def make_aberration_table(r,tol,write_csv,csv_file):
  """
  Determine a table of optical aberration angles for an observer in the Schwarzschild spacetime.
  Each line of the table is in the format [r,alpha,beta,beta-alpha,f], where r is observer's coordinate
  in units of Schwarzschild radius, alpha is angle as seen by observer if they're in a standard
  state of motion (free fall from rest at infinity), beta is angle on the celestial sphere,
  and f is the factor by which intensities are amplified by lensing. If this is the python version
  running, and write_csv is true, then we write a copy of the data to csv_file.
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
        if z<EPS and z>-10.0*EPS:
          z=EPS # tiny neg z happens due to rounding; make z big enough so it doesn't cause div by zero below
        if z<0.0:
          break # a photon with angular momentum this big can't exist at this r; won't happen with current algorithm
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
        PRINT("Geodesic at alpha=",alpha*180.0/MATH_PI," deg. is incomplete, done.")
        break
      w = info['user_data'] # winding number
      beta = atan2(final_x[3],final_x[2]) # returns an angle from -pi to pi
      if beta<0.0:
        beta = beta+2.0*MATH_PI
      beta = beta+w*2.0*MATH_PI
      table.append([r,alpha,beta])
      if verbosity>=2:
        PRINT("r=",r,", alpha=",alpha*180.0/MATH_PI," deg, beta=",beta*180.0/MATH_PI," deg")
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
    f = abs(sin(alpha)*dalpha/(sin(beta)*dbeta))
    if verbosity>=2:
      print("alpha=",alpha,", beta=",beta,", dalpha=",dalpha,", dbeta=",dbeta,", f=",f)
    # =dOmega(obs)/dOmega(infinity)=amplification, by Liouville's thm (?)
    # is abs() right?; beta can be >pi
    x.append(beta-alpha)
    x.append(f)
    table2.append(x)
#if "LANG" eq "python"
  if write_csv:
    with open(csv_file, 'w') as f:
      for x in table2:
        # x = [r,alpha,beta,beta-alpha,f]
        f.write(array_to_csv(x)+"\n")
    print("table of aberration data written to "+csv_file)
#endif
  return table2

def array_to_csv(x):
  return ",".join(map(lambda u : io_util.fl_n_decimals(u,12), x))

def read_csv_file(filename):
  with open(filename, 'r') as f:
    data = []
    reader = csv.reader(f)
    for row in reader:
      data.append(row)
    return data

def make_star_table(star_catalog,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,write_csv,csv_file):
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
  cursor = db.cursor()
  cursor.execute('''select max(id) from stars where mag<=?''',(max_mag,))
  n_stars = cursor.fetchone()[0]
  print("n_stars=",n_stars)
  m = euclidean.rotation_matrix_observer_to_celestial(ra_out,dec_out,1.0)
  m_inv = euclidean.rotation_matrix_observer_to_celestial(ra_out,dec_out,-1.0)
  table = []
  drawn = {}
  count_drawn = 0
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
      # Because the dec and RA windows are not quite rectangular, widen them a tiny bit to make sure
      # we don't miss any stars:
      k = 1.0+3.0/len(aberration_table)
      # more of a fudge factor if table is shorter; the 3 is empirical; we still miss 1 star out of 1600
      # with this fudge factor, and increasing the 3 to 20 doesn't fix that
      ddec = max_dec-min_dec
      dec = 0.5*(min_dec+max_dec)
      min_dec = dec-0.5*k*ddec
      max_dec = dec+0.5*k*ddec
      dra = max_ra-min_ra
      ra = 0.5*(min_ra+max_ra)
      min_ra = ra-0.5*k*dra
      max_ra = ra+0.5*k*dra
      # Estimate cut-off for brightness:
      f_approx = ab1[4]
      mag_corr = 2.5*log(f_approx)/log(10.0) # for f_approx<1, this is negative
      # Do a database search within this range of RA and dec.
      cursor = db.cursor()
      cursor.execute('''SELECT id,ra,dec,mag,bv FROM stars WHERE mag<=? AND ra>=? AND ra<=? AND dec>=? AND dec<=?''',\
                      (max_mag+mag_corr,min_ra,max_ra,min_dec,max_dec))
      all_rows = cursor.fetchall()
      for row in all_rows:
        id,ra,dec,mag,bv = row
        if id in drawn:
          continue
        is_rigel = abs(ra-rigel_ra_dec()[0])<0.01 and abs(dec-rigel_ra_dec()[1])<0.01 and mag<1.0
        if not is_rigel:
          continue # qwe
        drawn[id] = TRUE
        count_drawn = count_drawn+1
        beta,phi = celestial_to_beta(ra,dec,m_inv)
        alpha = alpha1+((alpha2-alpha1)/(beta2-beta1))*(beta-beta1)
        f = ab1[4]+((ab2[4]-ab1[4])/(beta2-beta1))*(beta-beta1) # interpolate amplification factor
        if f<EPS: # can have f<0 due to interpolation
          f=EPS
        brightness = exp(-(mag/2.5)*log(10.0)+log(f))
        table.append([alpha,phi,brightness,bv])
  db.close()
#if "LANG" eq "python"
  print("stars processed=",count_drawn," out of ",n_stars," with magnitudes under ",max_mag)
  if write_csv:
    with open(csv_file, 'w') as f:
      for x in table:
        # x = [alpha,phi,brightness,bv]
        f.write(array_to_csv(x)+"\n")
    print("table of star data written to "+csv_file)
#endif

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
  p2 = euclidean.apply_matrix(m,p)
  ncp = [0.0,0.0,1.0] # north celestial pole
  cx = [1.0,0.0,0.0] # x axis in celestial coords
  cy = [0.0,1.0,0.0] # x axis in celestial coords
  theta = acos(euclidean.dot(p2,ncp))
  phi = atan2(euclidean.dot(p2,cy),euclidean.dot(p2,cx))
  dec = 0.5*MATH_PI-theta
  ra = copy.copy(phi)
  if ra<0.0:
    ra = ra + 2.0*MATH_PI
  return [ra,dec]

def celestial_to_beta(ra,dec,m_inv):
  """
  Do the inverse of the transformation done by beta_to_celestial().
  """
  phi = copy.copy(ra)
  theta = 0.5*MATH_PI-dec
  p = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
  p2 = euclidean.apply_matrix(m_inv,p)
  zenith = [0.0,0.0,1.0]
  ox = [1.0,0.0,0.0] # x axis in obs. coords
  oy = [0.0,1.0,0.0] # x axis in obs. coords
  beta = 0.5*MATH_PI-acos(euclidean.dot(p2,zenith))
  phi = atan2(euclidean.dot(p2,oy),euclidean.dot(p2,ox))
  if beta<0.0:
    beta = -beta
    phi = phi+MATH_PI
  if phi>2.0*MATH_PI:
    phi = phi-2.0*MATH_PI
  if phi<0.0:
    phi = phi+2.0*MATH_PI
  return [beta,phi]

def rotation_matrix_observer_to_celestial(ra,dec,direction):
  # See beta_to_celestial() for explanation of what's going on.
  # direction = +-1; if direction is -1, compute the opposite transformation
  # Find the zenith vector, in celestial coords.
  zenith = [cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec)]
  # North celestial pole:
  ncp = [0.0,0.0,1.0]
  # Angle of rotation:
  rot = direction*acos(euclidean.dot(ncp,zenith))
  # Axis of rotation:
  if rot==0.0:
    axis = ncp # doesn't matter, just avoid division by zero that would otherwise happen in this case
  else:
    axis = euclidean.normalize(euclidean.cross_prod(zenith,ncp))
  return euclidean.rotation_matrix_from_axis_and_angle(rot,axis)

def rigel_ra_dec():
  """
  Returns the RA and dec of Rigel, in radians. Used for testing.
  """
  # 051432.27 -081205.9 +000001.9-000000.600004.2 00.18-0.03B8 0 0.05   2.07, bet Ori, Rigel
  # https://en.wikipedia.org/wiki/Rigel
  ra = ((5.0+14.0/60.0+32.27/3600.0)/24.0)*2*MATH_PI
  dec = (-(8+12/60.0+5.9/3600.)/360.0)*2*MATH_PI
  return [ra,dec]
  

main()
