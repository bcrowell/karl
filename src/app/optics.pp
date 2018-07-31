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

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild,euclidean,celestial,math_util,\
       star_properties,render
#if "LANG" eq "python"
import sys,os,copy,sqlite3,csv,random
import PIL,ephem
from PIL import Image
#endif

sloppy = TRUE
# ... I have checks at various places in the code to make sure that, e.g., a vector really is normalized
#     after I've constructed it in a way that should make it normalized. Sometimes these checks fail
#     due to rounding, e.g., near the horizon. Setting sloppy=TRUE disables these checks, so that
#     once the code is tested and I'm pretty sure it's correct, I no longer get exceptions being
#     thrown once in a while.


def main():
  if TRUE:
    r = 10.0
    if_fake = TRUE
    star_catalog_max_mag = 7
    width,height,fov_deg,view_rot_deg = [1200,600,130,100]
    do_image(r,"stars.png",if_fake,star_catalog_max_mag,width,height,fov_deg,view_rot_deg)
  else:
    # animation
    if_fake = FALSE # random stars wouldn't be at fixed locations
    star_catalog_max_mag = 9 # use dim stars to try to make up for lack of fake stars
    n = 100
    for i in range(n):
      z = float(i)/float(n)
      r = 9.0-7.9*z*z # accelerating, but not real physics for the state of motion we're using
      outfile = "animation"+("%03d" % i)+".png"
      print("---------------------- r=",r,", file=",outfile," ------------------------------")
      width,height,fov_deg,view_rot_deg = [600,300,90,100]
      do_image(r,outfile,if_fake,star_catalog_max_mag,width,height,fov_deg,view_rot_deg)

def do_image(r,image_file,if_fake,star_catalog_max_mag,width,height,fov_deg,view_rot_deg):
  verbosity=1
  star_catalog = '/usr/share/karl/mag'+str(star_catalog_max_mag)+'.sqlite'
  # Star catalog is built by a script in the directory data/star_catalog, see README in that directory.
  # falling inward from the direction of Rigel, https://en.wikipedia.org/wiki/Rigel :
  ra_out,dec_out = celestial.rigel_ra_dec()
  #ra_out,dec_out = celestial.antipodes_of_ra_and_dec(celestial.rigel_ra_dec())
  #ra_out,dec_out = celestial.antipodes_of_ra_and_dec(celestial.lmc_ra_dec())
  tol = 1.0e-3
  if_black_hole = TRUE
  max_mag = 12
  draw_sky(r,ra_out,dec_out,tol,"aberration.csv","v.csv","stars.csv","stars.json",verbosity,if_black_hole,\
            if_fake,max_mag,star_catalog,star_catalog_max_mag,width,height,fov_deg,view_rot_deg)
  os.system("src/render/render.rb stars.json "+image_file)

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
  aberration_table,v_table = make_aberration_tables(r,tol,verbosity)
  write_csv_file(aberration_table,aberration_csv,TRUE,"Table of aberration data written to")
  write_csv_file(v_table,v_table_csv,TRUE,"Table of ray velocities written to")
  star_table = make_star_table(star_catalog,aberration_table,v_table,r,if_black_hole,if_fake,\
                               ra_out,dec_out,max_mag,star_catalog_max_mag)
  write_csv_file(star_table,stars_csv,TRUE,"Table of star data written to")
  render.render(star_table,image_json,verbosity,width,height,fov_deg,view_rot_deg)

#--------------------------------------------------------------------------------------------------

def make_aberration_tables(r,tol,verbosity):
  """
  Determine a table of optical aberration angles for an observer in the Schwarzschild spacetime.
  Each line of the table is in the format [r,alpha,beta,beta-alpha,f].
  Here r is observer's coordinate
  in units of Schwarzschild radius, alpha is angle as seen by observer if they're in a standard
  state of motion (free fall from rest at infinity), beta is angle on the celestial sphere,
  and f is the factor by which intensities are amplified by lensing. 
  In addition, we return a separate table [alpha,v0,...v4] of the 
  components of the velocity vector v^kappa of the light ray at observation,
  in 5-dimensional Schwarzschild coordinates. 
  This is used in order to find the Doppler shifts later.
  The alpha values in this table are a subset of the
  ones in the aberration table, because it seems safer to interpolate the Doppler shifts at the end,
  rather than interpolating the vector components here.
  The normalization of the velocities is such that the velocity vector at emission is (1,-1,0,0,0).
  The observer's velocity is not calculated here because it can easily be found later, and
  would be the same for all rays: -- u_kappa=(1,-A*sqrt(1-A)), where A=1-1/r. Note that this
  is the covariant form of u.
  """
  is_schwarzschild = TRUE
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  x_obs,v_obs,rho = schwarzschild_standard_observer(r,spacetime,chart,pars)
  if r>1.0:
    max_le = r/sqrt(aa) # maximum possible L/E for a photon at this r (not the max. for visibility)
  else:
    max_le = 10.0 # fixme
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
  v_table = []
  last_deflection = 0.0
  s0 = 2
  # ...scaling factor for number of angular steps
  #    Making this a big value, like 10, makes fake stars take a long time, because the sky is subdivided
  #    very finely. Making it a very small value, like 2, may make interpolation of Doppler shifts too crude.
  s1 = 2 # for moderate deflections, reduce step size by this factor
  s2 = 10 # ... additional factor for large deflections
  n_angles = s0*s1*s2 # we actually skip most of these angles
  for in_n_out in range(2): # 0=outward, 1=inward
    for i in range(n_angles):
      frac_skip = s1*s2 # normally we fill in most of the angles by interpolation
      if last_deflection>0.5:
        frac_skip = frac_skip//s1
      if last_deflection>2.0:
        frac_skip = frac_skip//s2
      skip_this = not (i%frac_skip==0)
      z = (float(i)/float(n_angles)) # note that this should *not* be n_angles-1
      if in_n_out==0:
        le = z*max_le
      else:
        if is_schwarzschild:
          le_ph = 0.5*3**1.5 # L/E of the photon sphere, unstable circular orbits for photons in Schwarzschild
          le = math_util.linear_interp(0.0,1.0,max_le,le_ph,z)
        else:
          # more generic logic, more likely to be appropriate if not Sch.
          le = (1.0-z)*max_le
      alpha,v_observation = le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
      x = x_obs # initial position
      v = CLONE_ARRAY_OF_FLOATS(v_observation)
      # ... Initial velocity of ray, which is simulated going backward in time, as if emitted rather than absorbed.
      #     Will get rescaled later, see comments below, but clone it to make sure it can't get munged.
      if skip_this:
        beta,done = [NONE,FALSE] # fill in later by interpolation
      else:
        beta,done,v_emission = do_ray(spacetime,chart,pars,x,v,r,tol,count_winding,alpha)
      table.append([r,alpha,beta])
      got_result = not (IS_NONE(beta) or IS_NAN(beta))
      #---
      # Calculate the velocity of the ray at observation. This would actually be pretty trivial, since
      # v is the "initial" velocity for solving the diffeq, but is actually the final
      # velocity of the ray, at observation. However, we do need v_emission for normalization.
      # Because the Schwarzschild coordinates are spherical coordinates, v_emission is guaranteed
      # to have the form v^kappa=(q,-q,0,0,0). We need to retrofit v=v_observation to have a normalization
      # that would have resulted from emission with v^kappa=(1,-1,0,0,0), which corresponds to a different  
      # choice of affine parameter. In the Schwarzschild case, v_t is conserved, so the result for
      # v_observation should always be v^t=1/(1-1/r), and we would not need the results of ray-tracing
      # to tell us this, but I want this code to be more general, so I don't use that assumption.
      if got_result:
        q = v_emission[0]
        v_observation = vector.scalar_mult(v_observation,1.0/q)
        # This vector was set up for an outgoing ray, but we're flipping this around, treating it as an incoming
        # ray that was observed here. Therefore we need to flip components 1...4:
        v_observation = vector.scalar_mult(v_observation,-1.0)
        v_observation[0] = -v_observation[0] # flip 0 component back to what it was
        v_table.append([alpha,]+v_observation)
        # ... The + here is concatenation of the lists. This won't work in js.
      if got_result:
        if verbosity>=2:
          if alpha==0.0:
            err_approx = 0.0
          else:
            approx = alpha/(sqrt(r)-1)
            err_approx = (approx-abs(alpha-beta))/alpha
            PRINT("r=",r,", alpha=",alpha*180.0/MATH_PI," deg, beta=",beta*180.0/MATH_PI,\
                " deg, rel err in approx=",err_approx)
        last_deflection = abs(alpha-beta)
        if abs(alpha-beta)>5.0*MATH_PI: # Riazuelo says 5pi is enough to get all visual effects.
          PRINT("Deflection=",abs(alpha-beta)*180.0/MATH_PI," deg. is >5pi, done.")
          done = TRUE
      if done:
        break
    if done:
      break
  fill_in_aberration_table_by_interpolation(table,r)
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
  return [table2,v_table]

def do_ray(spacetime,chart,pars,x,v,r,tol,count_winding,alpha):
  n = 100
  ndebug=0
  if verbosity>=3:
    ndebug=n/10
  ri = r
  lambda_max = 10.0 # fixme, sort of random
  info = {}
  while True:
    opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,'future_oriented':FALSE,'tol':tol,
          'user_function':count_winding,'no_enter_horizon':TRUE}
    err,final_x,final_v,final_a,final_lambda,info,sigma  = \
            fancy.trajectory_schwarzschild(spacetime,chart,pars,x,v,opt)
    if err!=0:
      if err==RK_INCOMPLETE or err==RK_TRIGGER:
        PRINT("Geodesic at alpha=",alpha*180.0/MATH_PI," deg. is incomplete or entered horizon, done.")
        return [NAN,TRUE,NONE]
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
  w = info['user_data'] # winding number
  beta = atan2(final_x[3],final_x[2]) # returns an angle from -pi to pi
  if beta<0.0:
    beta = beta+2.0*MATH_PI
  beta = beta+w*2.0*MATH_PI
  return [beta,FALSE,final_v]

def fill_in_aberration_table_by_interpolation(table,r):
  #-------------------
  # Use interpolation to fill in missing values.
  # Find i and j that both have real data.
  for i in range(len(table)-2):
    if IS_NONE(table[i][2]):
      continue
    j=i+1
    while j<=len(table)-1 and IS_NONE(table[j][2]):
      j=j+1
    if j==i+1 or j>len(table)-1:
      continue
    alpha1 = table[i][1]
    alpha2 = table[j][1]
    beta1 = table[i][2]
    beta2 = table[j][2]
    # Find error compared to small-angle approx:
    err1 = alpha1/(sqrt(r)-1)-(beta1-alpha1)
    err2 = alpha2/(sqrt(r)-1)-(beta2-alpha2)
    k=i+1
    while IS_NONE(table[k][2]) and k<j:
      alpha = table[k][1]
      err = math_util.linear_interp(alpha1,alpha2,err1,err2,alpha) # interpolate to find est. of what error would have been
      beta = alpha+alpha/(sqrt(r)-1)-err
      table[k][2] = beta
      k=k+1
  # Remove any uncomputed items from the end of the table:
  while IS_NONE(table[len(table)-1][2]):
    table.pop()
  while IS_NAN(table[len(table)-1][2]):
    table.pop()

def schwarzschild_standard_observer(r,spacetime,chart,pars):
  """
  The main input is the Schwarzschild radial coordinate r of an observer in the Schwarzschild spacetime.
  The other inputs are described in comments in the calling code.
  Although the inputs include data about the spacetime and chart, this code will not actually work except
  for Schwarzschild cooordinates in the Schwarzschild spacetime.
  The output is the three vectors [x_obs,v_obs,rho] describing an observer at the standard state of motion,
  free-falling from rest at infinity. These are the observer's position and velocity, and
  the observer's radial vector rho, parallel-transported in from infinity.
  """
  aa = 1-1/r
  x_obs = [0.0,r,1.0,0.0,0.0]
  # Calculate the velocity vector of the observer, for the standard state of motion:
  v_obs = [1/aa,-sqrt(1-aa),0.0,0.0,0.0]
  #    ... solution of the equations E=Atdot=1, |v|^2=Atdot^2-(1/A)rdot^2=1
  #        We check below that it actually is a solution.
  # Calculate the observer's radial vector rho, parallel-transported in from infinity.
  rho = vector.proj(spacetime,chart,pars,x_obs,v_obs,[0.0,1.0,0.0,0.0,0.0])
  # ... take rhat and project out the part orthogonal to the observer's velocity; not yet normalized
  rho = vector.normalize_spacelike(spacetime,chart,pars,x_obs,rho)
  # ... now normalize it to have |rho|^2=-1
  if abs(vector.norm(spacetime,chart,pars,x_obs,rho)+1.0)>EPS*10 and not sloppy:
    THROW("norm(rho)!=-1")
  energy_obs = aa*v_obs[0]
  if abs(energy_obs-1.0)>EPS*10:
    THROW("energy of observer="+str(energy_obs))
  if abs(vector.norm(spacetime,chart,pars,x_obs,v_obs)-1.0)>EPS*10 and not sloppy:
    THROW("observer's velocity has bad norm")
  return [x_obs,v_obs,rho]

def alpha_max_schwarzschild(r):
  """
  Maximum alpha at which an observer at r, in the standard state of motion, receives rays.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  le_ph = 0.5*3**1.5 # L/E of the photon sphere, unstable circular orbits for photons in Schwarzschild
  x_obs,v_obs,rho = schwarzschild_standard_observer(r,spacetime,chart,pars)
  # In the following, there is no reverse ray-tracing going on, so inward means the observer really
  # absorbs the photon at a point where it has dr/dt<0.
  # in_n_out is 0 if outward, 1 if inward
  if r>=1.5:
    in_n_out = 0 # ray peeled off of the photon sphere and came outward to us
  else:
    in_n_out = 1 # limiting rays peel off of the photon sphere and fall inward to us
  alpha_max,v_observation = le_to_alpha_schwarzschild(r,le_ph,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars)
  return alpha_max

def le_to_alpha_schwarzschild(r,le,in_n_out,x_obs,v_obs,rho,spacetime,chart,pars):
  """
  The main inputs are Schwarzschild radius r of the observer, the ratio L/E of the angular momentum of a photon
  to its energy, and in_n_out, which equals 0 if the photon is going radially outward, 1 if inward.
  The other inputs are described in comments in the calling code.
  Although the inputs include data about the spacetime and chart, this code will not actually work except
  for Schwarzschild cooordinates in the Schwarzschild spacetime.
  Output is [alpha,v_observation], where v_observation is the velocity vector of the observer
  in the standard state of motion (infalling from rest at infinity), and alpha is the azimuthal angle of
  the photon as measured by the observer.
  """
  aa = 1-1/r
  #----
  # Find velocity vector of the photon.
  if le==0.0:
    v = [1.0,aa,0.0,0.0,0.0] # should probably change sign of v_r based on in_n_out, but test carefully -- qwe
  else:
    z = r*r/(le*le)-aa
    if z<EPS and z>-10.0*EPS:
      z=EPS # tiny neg z happens due to rounding; make z big enough so it doesn't cause div by zero below
    if z<0.0:
      THROW("z<0, a photon with angular momentum this big can't exist at this r")
      # ... won't happen with current algorithm
    dphi_dr = (1.0/r)/sqrt(z)
    if in_n_out==1:
      dphi_dr = -dphi_dr
    dphi_dt = le*aa/(r*r)
    v = [1.0,dphi_dt/dphi_dr,0.0,dphi_dt,0.0] # tangent vector
  #----
  # Check that its norm is zero.
  norm = vector.norm(spacetime,chart,pars,x_obs,v)
  if abs(norm)>EPS*10 and not sloppy:
    THROW("norm="+str(norm))
  #----
  v_perp = vector.proj(spacetime,chart,pars,x_obs,v_obs,v)
  # ... part of photon's velocity orthogonal to observer's velocity
  v_perp = vector.normalize_spacelike(spacetime,chart,pars,x_obs,v_perp)
  # ... normalized
  zzz = math_util.force_into_range(-vector.inner_product(spacetime,chart,pars,x_obs,rho,v_perp),-1.0,1.0)
  # ... Minus sign is because this is +--- signature, so euclidean dot products have flipped signs.
  #     I thought I convinced myself that there should be another - because we're doing direction
  #     ray appears to have come from, not direction it's going, but that seems to mess things up.
  alpha = acos(zzz)
  # ... angle at which the observer says the photon is emitted, see docs; the force_into_range() is
  #     necessary because sometimes we get values for zzz like 1.0000000000000002 due to rounding
  return [alpha,v]

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
  The star table returned by this routine is in the format [].
  """
  m = celestial.rotation_matrix_observer_to_celestial(ra_out,dec_out,1.0)
  m_inv = celestial.rotation_matrix_observer_to_celestial(ra_out,dec_out,-1.0)
  table = []
  table,stats_real =\
           real_stars(table,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,star_catalog,m,m_inv,\
                      v_table)
  n_stars = stats_real['n_stars']
  count_drawn = stats_real['count_drawn']
  PRINT("stars processed=",count_drawn," out of ",n_stars," with apparent magnitudes under ",max_mag)
  if if_fake:
    table,stats_fake =\
           fake_stars(table,aberration_table,r,if_black_hole,ra_out,dec_out,max_mag,m,m_inv,star_catalog_max_mag,\
                      v_table)
    count_fake = stats_fake['count_fake']
    PRINT("Drew ",count_fake," fake stars.")
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
    alpha = math_util.linear_interp_from_table(aberration_table,2,1,beta,0,len(aberration_table)-1)
    count_drawn = count_drawn+1
    if if_black_hole:
      f = math_util.linear_interp_from_table(aberration_table,2,4,beta,0,len(aberration_table)-1)
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
            beta = math_util.linear_interp(alpha1,alpha2,beta1,beta2,alpha)
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
  v0 = math_util.linear_interp_from_table(v_table,0,1,alpha,0,len(v_table)-1)
  v1 = math_util.linear_interp_from_table(v_table,0,2,alpha,0,len(v_table)-1)
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
