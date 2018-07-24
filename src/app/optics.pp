#!/usr/bin/python3
#####################

# This is only going to work in python, not javascript.
# Requires debian packages python3-pil python3-ephem

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild,euclidean,celestial
#if "LANG" eq "python"
import sys,os,copy,sqlite3,csv,random
import PIL,ephem
from PIL import Image
#endif

verbosity=2

star_catalog_max_mag = 7
star_catalog = '/usr/share/karl/mag7.sqlite'
# Star catalog is built by a script in the directory data/star_catalog, see README in that directory.

def main():
  r = 9.0
  # falling inward from the direction of Rigel, https://en.wikipedia.org/wiki/Rigel :
  ra_out,dec_out = celestial.rigel_ra_dec()
  tol = 1.0e-3
  aberration_table = make_aberration_table(r,tol,TRUE,"aberration.csv")
  max_mag = 12 # max apparent mag to show; causes random fake stars to be displayed down to this mag,
              # in addition to the ones bright enough to be in the catalog
  star_table = make_star_table(star_catalog,aberration_table,r,TRUE,ra_out,dec_out,max_mag,TRUE,"stars.csv")

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
  last_deflection = 0.0
  s1 = 5 # for moderate deflections, reduce step size by this factor
  s2 = 20 # ... additional factor for large deflections
  n_angles = 10*s1*s2 # we actually skip most of these angles
  for in_n_out in range(2): # 0=outward, 1=inward
    for i in range(n_angles):
      frac_skip = s1*s2 # normally we skip most of of the angles, don't need the resolution
      if last_deflection>0.5:
        frac_skip = frac_skip//s1
      if last_deflection>2.0:
        frac_skip = frac_skip//s2
      if not (i%frac_skip==0):
        continue
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
        # need more precision for rays that are greatly deflected
        tt = tol
        if last_deflection>0.3:
          tt = tt*0.001 
        if last_deflection>1.0:
          tt = tt*0.001 
        opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,'tol':tt,
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
      last_deflection = abs(alpha-beta)
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
  m = celestial.rotation_matrix_observer_to_celestial(ra_out,dec_out,1.0)
  m_inv = celestial.rotation_matrix_observer_to_celestial(ra_out,dec_out,-1.0)
  table = []
  drawn = {}
  count_drawn = 0
  count_fake = 0
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
    dbeta = beta2-beta1
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
          ra,dec = celestial.beta_to_celestial(beta_c,phi_c,m)
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
      count_found = 0
      for row in all_rows:
        id,ra,dec,mag,bv = row
        if id in drawn:
          continue
        count_found = count_found+1
        drawn[id] = TRUE
        count_drawn = count_drawn+1
        beta,phi = celestial.celestial_to_beta(ra,dec,m_inv)
        alpha = alpha1+((alpha2-alpha1)/(beta2-beta1))*(beta-beta1)
        brightness = brightness_helper(beta,mag,ab1[4],ab2[4],beta1,beta2)
        table.append([alpha,phi,brightness,bv])
      if star_catalog_max_mag<max_mag+mag_corr:
        # Fill in fake random stars too dim to have been in the catalog.
        # Seares, 1925, http://adsbit.harvard.edu//full/1925ApJ....62..320S/0000320.000.html
        # Frequency of magnitude m is propto exp(am), where a=0.86 (my estimate from their graphs).
        # This means that if the number of stars up to magnitude m0 is N, then the number from
        # m0 to m0+1 is (e^a-1)N=1.36N.
        # Could also use GAIA full sky maps, either to estimate the density of stars or
        # to provide a pretty milky way background, but the latter would not be something we could
        # easily simulate doppler shifts of.
        #   http://sci.esa.int/gaia/60196-gaia-s-sky-in-colour-equirectangular-projection/
        #   http://sci.esa.int/gaia/60170-gaia-s-new-map-of-star-density/
        nn = 27.1
        # ... normalization of GAIA density map, found empirically by matching to results from
        #     star atlas
        solid_angle = abs(dphi*sin(beta)*dbeta)
        #n = count_found
        n = get_star_density(ra,dec)*nn*solid_angle
        k=CEIL(max_mag+mag_corr-star_catalog_max_mag) # number of missing magnitudes to simulate
        for i in range(k):
          dn = 1.36*n
          n = n+dn
          for j in range(poisson_random(dn)):
            count_fake = count_fake+1
            alpha = uniform_random(alpha1,alpha2)
            phi = uniform_random(phi1,phi2)
            mag = uniform_random(star_catalog_max_mag+i,star_catalog_max_mag+i+1)
            beta = beta1+((beta2-beta1)/(alpha2-alpha1))*(alpha-alpha1)
            brightness = brightness_helper(beta,mag,ab1[4],ab2[4],beta1,beta2)
            bv = 0.0 # fixme
            if brightness>1.0e-3 and alpha <1.5:
              print("very bright fake star at low alpha, alpha,phi,brightness,mag=",alpha,phi,brightness,mag)
            table.append([alpha,phi,brightness,bv])
  db.close()
#if "LANG" eq "python"
  print("stars processed=",count_drawn," out of ",n_stars," with apparent magnitudes under ",max_mag,\
            ", plus ",count_fake," fake stars")
  if write_csv:
    with open(csv_file, 'w') as f:
      for x in table:
        # x = [alpha,phi,brightness,bv]
        f.write(array_to_csv(x)+"\n")
    print("table of star data written to "+csv_file)
#endif

def uniform_random(a,b):
  return a+(b-a)*random.random()

def poisson_random(mean):
  # rough approximation to a Poisson random variable
  if mean<0.3:
    if random.random()<mean:
      return 1
    else:
      return 0
  if mean<3.0:
    # https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
    l = exp(-mean)
    k=0
    p=1
    while TRUE:
      k = k+1
      u = random.random()
      p = p*u
      if p<l:
        break
    return k-1
  # large mean
  x = int(random.gauss(mean,sqrt(mean))+0.5)
  if x<0:
    x=0
  return x

def brightness_helper(beta,mag,f1,f2,beta1,beta2):
  f = f1+((f2-f1)/(beta2-beta1))*(beta-beta1) # interpolate amplification factor
  if f<EPS: # can have f<0 due to interpolation
    f=EPS
  brightness = exp(-(mag/2.5)*log(10.0)+log(f))
  return brightness

def get_star_density(ra,dec):
  # Normalized to 1 for densest part of map.
  # This is for use with the GAIA map, which is downloaded and chopped up by scripts in data/star_density.
  lon,lat = to_galactic_lon_lat(ra,dec)
  x_res = 4000 # pixels
  y_res = 2000
  x = CEIL((lon/(2*MATH_PI))*x_res)+x_res//2 # I assume they put zero galactic lon in the center.
  y = CEIL((lat/(MATH_PI/2.0))*y_res)+y_res//2
  if x>x_res-1:
    x = x-x_res
  x = force_into_range(x,0,x_res-1)
  y = force_into_range(y,0,y_res-1)
  # Filenames are like gaia_color_equirect_medium_part_08_07.png
  # Here, i=8, j=7. See data/star_density/split_gaia.rb.
  # Blocks are 100x100 pixels. The files are in equirectangular projection, https://en.wikipedia.org/wiki/Equirectangular_projection ,
  # which basically means that the lat-lon conversion is trivial.
  ww = 200
  hh = 200
  i = x//ww
  j = y//hh
  ii = "%02d" % i
  jj = "%02d" % j
  png = '/usr/share/karl'+'/gaia_color_equirect_medium_part_'+ii+'_'+jj+'.jpg'
  image = Image.open(png)
  pixels = image.load()
  r,g,b = pixels[x%ww,y%hh]
  image.close()
  result = float(r+g+b)/(255.0*3.0)
  return result

def to_galactic_lon_lat(ra,dec):
  """
  Inputs are RA and declination in equatorial coordinates, outputs are galactic lon and lat.
  Inputs and outputs are in radians.
  """
  g = ephem.Galactic(ephem.Equatorial(ra, dec))
  return [g.lon+0.0,g.lat+0.0]

def force_into_range(x,a,b):
  if x<a:
    return a
  if x>b:
    return b
  return x

main()
