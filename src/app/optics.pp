#!/usr/bin/python3
#####################

# This is only going to work in python, not javascript.
# Code that could also work in js should be moved into libraries such as star_properties and ray.
# Requires debian packages python3-pil python3-ephem

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild,euclidean,celestial,math_util,\
       star_properties,render,ray,animation
#if "LANG" eq "python"
import sys,os,copy,sqlite3,csv,random
import PIL,ephem,datetime
from PIL import Image
#endif


def main():
  do_what = 3
  if do_what==1:
    r = 0.9
    alpha = 2.3754638 # very close to alpha_max=2.3759564949418355
    tol = 1.0e-3
    alpha_max = ray.alpha_max_schwarzschild(r)
    print("r=",r,", alpha_max=",alpha_max)
    count_winding(0.0,[],[],0,0,{})
    beta,if_incomplete,final_v = ray.do_ray_schwarzschild(r,tol,count_winding,alpha)
    print("alpha=",alpha,", beta=",beta)
    exit(0)
  if do_what==2:
    # Make one image.
    r = 0.9
    if_fake = TRUE
    star_catalog_max_mag = 7
    width,height,fov_deg,view_rot_deg = [1200,600,130,100]
    verbosity = 3
    do_image(r,"stars.png",if_fake,star_catalog_max_mag,width,height,fov_deg,view_rot_deg,"","")
    exit(0)
  if do_what==3:
    # animation
    print("Control the following parameters  in makefile with command-line args:")
    prep_level = int(sys.argv[1])
    segment = int(sys.argv[2])
    if len(sys.argv)>=4:
      single_frame = int(sys.argv[3])
    else:
      single_frame = -1
    print("  prep level=",prep_level)
    print("  segment=",segment)
    print("  single frame=",single_frame)
    # prep levels:
    #   1 -- no images generated, just helps with planning motion of observer
    #   2 -- initial and final frames only, stars only down to magnitude 7
    #   3 -- initial, final, and every 10th frame
    #   4 -- real thing
    animation.do_frames(segment,prep_level,single_frame,do_image)

def count_winding(lam,x,v,spacetime,chart,pars):
  # not normally needed, just used by some test code
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

def do_image(r,image_file,if_fake,star_catalog_max_mag,width,height,fov_deg,view_rot_deg,\
                        data_file_prefix,data_file_suffix):
  verbosity=3
  star_catalog = '/usr/share/karl/mag'+str(star_catalog_max_mag)+'.sqlite'
  # Star catalog is built by a script in the directory data/star_catalog, see README in that directory.
  # falling inward from the direction of Rigel, https://en.wikipedia.org/wiki/Rigel :
  ra_out,dec_out = celestial.rigel_ra_dec()
  #ra_out,dec_out = celestial.antipodes_of_ra_and_dec(celestial.rigel_ra_dec())
  #ra_out,dec_out = celestial.antipodes_of_ra_and_dec(celestial.lmc_ra_dec())
  tol = 1.0e-3
  if_black_hole = TRUE
  max_mag = 12
  aberration_csv = data_file_prefix+"aberration"+data_file_suffix+".csv"
  v_csv          = data_file_prefix+"v"+data_file_suffix+".csv"
  stars_csv      = data_file_prefix+"stars"+data_file_suffix+".csv"
  stars_json     = data_file_prefix+"stars"+data_file_suffix+".json"
  draw_sky(r,ra_out,dec_out,tol,aberration_csv,v_csv,stars_csv,stars_json,verbosity,if_black_hole,\
            if_fake,max_mag,star_catalog,star_catalog_max_mag,width,height,fov_deg,view_rot_deg)
  os.system("src/render/render.rb "+stars_json+" "+image_file)

#--------------------------------------------------------------------------------------------------

def draw_sky(r,ra_out,dec_out,tol,aberration_csv,v_table_csv,stars_csv,image_json,verbosity,\
             if_black_hole,if_fake,max_mag,star_catalog,star_catalog_max_mag,width,height,\
             fov_deg,view_rot_deg):
  """
  Draw the sky as seen by an observer near a Schwarzschild black hole.
  The observer is at radius r and angular coordinates specified by the given
  RA and dec, i.e., the black hole is placed at the center of the earth's celestial sphere.
  Use tolerance tol for geodesics. The two csv files are currently only for the user to inspect,
  if desired. The output is a set of pixel arrays written to image_json, which is actually
  rendered to a png file by render.rb.
  The variable max_mag is the max apparent mag to show; causes random fake stars to be displayed down
  to this mag, in addition to the ones bright enough to be in the catalog.
  If if_black_hole is false, the star field is drawn as it would appear without any black hole.
  Width and height are in pixels.
  """
  max_deflection = 5.0*MATH_PI  # ... Riazuelo says 5pi is enough to get all important effects.
  aberration_table,v_table = ray.make_aberration_tables(r,tol,verbosity,max_deflection)
  sys.stdout.flush()
  write_csv_file(aberration_table,aberration_csv,TRUE,"Table of aberration data written to")
  write_csv_file(v_table,v_table_csv,TRUE,"Table of ray velocities written to")
  star_table = make_star_table(star_catalog,aberration_table,v_table,r,if_black_hole,if_fake,\
                               ra_out,dec_out,max_mag,star_catalog_max_mag)
  sys.stdout.flush()
  write_csv_file(star_table,stars_csv,TRUE,"Table of star data written to")
  render.render(star_table,image_json,verbosity,width,height,fov_deg,view_rot_deg)
  sys.stdout.flush()

def make_star_table(star_catalog,aberration_table,v_table,r,if_black_hole,if_fake,ra_out,dec_out,max_mag,\
                    star_catalog_max_mag):
  """
  Make a table of stars seen by on observer in the vicinity of a black hole, in a standard state of motion,
  which is free fall from rest at infinity.
  If if_black_hole is false, then the black hole is omitted.
  The parameters ra_out,dec_out are the RA and declination, in radians, of the point on the celestial
  sphere from which the observer has fallen.
  The aberration_table is generated by make_aberration_table(), and contains rows in the format
  [r,alpha,beta,beta-alpha,f]. Here alpha is the observer's angle.
  Max_mag gives the maximum apparent magnitude to include.
  The star table returned by this routine is in the format [alpha,phi,brightness,log temp] (see
  star_table_entry_helper).
  """
  PRINT("making star table")
  sys.stdout.flush()
  m = celestial.rotation_matrix_observer_to_celestial(ra_out,dec_out,1.0)
  m_inv = celestial.rotation_matrix_observer_to_celestial(ra_out,dec_out,-1.0)
  table = []
  table,stats_real =\
           real_stars(table,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,star_catalog,m,m_inv,\
                      v_table)
  n_stars = stats_real['n_stars']
  count_drawn = stats_real['count_drawn']
  PRINT("  stars processed=",count_drawn," out of ",n_stars," with apparent magnitudes under ",max_mag)
  sys.stdout.flush()
  if if_fake:
    table,stats_fake =\
           fake_stars(table,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,m,m_inv,star_catalog_max_mag,\
                      v_table)
    count_fake = stats_fake['count_fake']
    PRINT("  Drew ",count_fake," fake stars.")
  PRINT("  Done making star table.")
  sys.stdout.flush()
  return table

def real_stars(table,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,star_catalog,m,m_inv,v_table):
  db = sqlite3.connect(star_catalog)
  cursor = db.cursor()
  cursor.execute('''select max(id) from stars''')
  n_stars = cursor.fetchone()[0]
  print("n_stars=",n_stars)
  count_drawn = 0
  for i in range(n_stars):
    cursor = db.cursor()
    cursor.execute('''SELECT ra,dec,mag,bv FROM stars WHERE id=?''',(i+1,))
    row = cursor.fetchall()[0]
    ra,dec,mag,bv = row
    if mag>max_mag:
      continue
    beta,phi = celestial.celestial_to_beta(ra,dec,m_inv)
    alpha,unsafe = math_util.linear_interp_from_table_safe(aberration_table,2,1,beta,FALSE,FALSE,0.0)
    if unsafe:
      PRINT("beta=",beta)
      THROW('unsafe extrapolation, this can happen if the grid for tabulating beta(alpha) is too coarse')
    count_drawn = count_drawn+1
    if if_black_hole:
      f = math_util.linear_interp_from_table(aberration_table,2,4,beta)
      if f<EPS: # can have f<0 due to interpolation
        f=EPS
    else:
      f=1.0
    brightness = exp(-(mag/2.5)*log(10.0)+log(f))
    ln_temp = star_properties.bv_to_log_temperature(bv)
    star_table_entry_helper(table,alpha,phi,brightness,ln_temp,beta,if_black_hole,r,v_table)
  db.close()
  PRINT("done with real stars, drew ",count_drawn," with magnitudes less than ",max_mag)
  return [table,{'n_stars':n_stars,'count_drawn':count_drawn}]

def fake_stars(table,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,m,m_inv,star_catalog_max_mag,\
                 v_table):
  count_fake = 0
  count_boxes = 0
  for i in range(len(aberration_table)-1):
    ab1 = aberration_table[i]
    ab2 = aberration_table[i+1]
    alpha1 = ab1[1]
    alpha2 = ab2[1]
    alpha = 0.5*(alpha1+alpha2) # midpoint of strip
    #PRINT("adding fake stars for alpha=",alpha*180.0/MATH_PI," deg., count_boxes=",count_boxes)
    dalpha = alpha2-alpha1
    beta1 = ab1[2]
    beta2 = ab2[2]
    beta = 0.5*(beta1+beta2) # midpoint
    dbeta = beta2-beta1
    # Break up the 2pi range of azimuthal angles into n_az chunks, so that each chunk within this
    # strip is approximately square.
    n_az = CEIL(2.0*MATH_PI*sin(alpha)/dalpha)
    count_boxes = count_boxes+n_az
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
            beta,unsafe = math_util.linear_interp_safe(alpha1,alpha2,beta1,beta2,alpha,FALSE,FALSE,0.0)
            if unsafe:
              THROW('unsafe extrapolation in aberration table')
            brightness = brightness_helper(beta,mag,ab1[4],ab2[4],beta1,beta2,if_black_hole)
            bv = star_properties.spectral_class_to_bv(star_properties.random_spectral_class())
            ln_temp = star_properties.bv_to_log_temperature(bv)
            star_table_entry_helper(table,alpha,phi,brightness,ln_temp,beta,if_black_hole,r,v_table)
  return [table,{'count_fake':count_fake}]

def star_table_entry_helper(table,alpha,phi,brightness,ln_temp,beta,if_black_hole,r,v_table):
  if if_black_hole:
    angle = alpha
    doppler = doppler_helper(r,v_table,alpha)
  else:
    angle = beta
    doppler = 1.0
  ln_temp_doppler_corr = log(doppler)
  # ... https://en.wikipedia.org/wiki/Black-body_radiation#Doppler_effect_for_a_moving_black_body
  #     A Doppler-shifted blackbody spectrum is still a blackbody spectrum, with a shifted temp.
  #     Doppler effect on whole-spectrum intensity is not handled here, was included in the
  #     factor f that we got from Liouville's theorem.
  #     Doppler shifts can shift energy into or out of the visible-light spectrum, but we don't do anything
  #     about that here. This is just physics, not vision.
  table.append([angle,phi,brightness,ln_temp+ln_temp_doppler_corr])

def doppler_helper(r,v_table,alpha):
  # omega'/omega = (u'_k v'^k)/(u_m v^m), where u=emitter or observer, v=ray, unprimed=emission, primed=observed.
  # The following contains some logic, clearly marked, that assumes Schwarzschild spacetime.
  # Ditto for assumptions about the observer's state of motion.
  denominator = 1.0
  # ... Denominator of the expression above for omega'/omega.
  #     Because Schwarzschild coordinates are spherical and asymptotically flat, we always have
  #     u_m=(1,0,zeroes) and v^m=(1,-1,zeroes), so this inner product is 1.
  aa = 1.0-1.0/r # assumes Sch.
  u0 = 1.0
  u1 = (1.0/aa)*sqrt(1-aa)
  # ... Covariant form of (t,r) components of observer's velocity. Assumes Sch. and observer radially infalling
  #     from rest at infinity. u0=1 because it's the conserved energy in Sch. spacetime and we started from
  #     rest at infinity. u1 follows from normalization.
  v0 = math_util.linear_interp_from_table(v_table,0,1,alpha)
  v1 = math_util.linear_interp_from_table(v_table,0,2,alpha)
  numerator = u0*v0+u1*v1 # inner product of u' (covariant) with v' (contravariant).
  doppler = numerator/denominator
  return doppler

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

def brightness_helper(beta,mag,f1,f2,beta1,beta2,if_black_hole):
  if if_black_hole:
    f = math_util.linear_interp(beta1,beta2,f1,f2,beta) # interpolate amplification factor
    if f<EPS: # can have f<0 due to interpolation
      f=EPS
  else:
    f=1.0
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
  x = math_util.force_into_range(x,0,x_res-1)
  y = math_util.force_into_range(y,0,y_res-1)
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


def array_to_csv(x):
  return ",".join(map(lambda u : io_util.fl_n_decimals(u,12), x))

def read_csv_file(filename):
#if "LANG" eq "python"
  with open(filename, 'r') as f:
    data = []
    reader = csv.reader(f)
    for row in reader:
      data.append(row)
    return data
#else
  return
#endif

def write_csv_file(table,filename,if_message,message):
#if "LANG" eq "python"
  with open(filename, 'w') as f:
    for x in table:
      f.write(array_to_csv(x)+"\n")
  if if_message:
    print(message+" "+filename)
#else
  return
#endif

main()
