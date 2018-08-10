#include "language.h"
#include "math.h"
#include "spacetimes.h"
#include "precision.h"

import sys,os,csv,datetime
import ray,render

def do_frames(segment,prep_level,do_image):
  width,height,fov_deg,view_rot_deg = [1280,720,110,80]
  # ... 1280x720 with square pixels is 16:9, recommended by youtube
  if_fake = FALSE # random stars wouldn't be at fixed locations
  if prep_level>=3:
    star_catalog_max_mag = 9 # use dim stars to try to make up for lack of fake stars
  else:
    star_catalog_max_mag = 7
  if segment==1:
    # Infalling from r=30 to 9, not super exciting, no need for high time resolution.
    n = 200
    r = 30.0
    dtau = 0.4607 # lands us at very close to r=9 on final frame
  if segment==2:
    # r=9 to 3, slightly higher time resolution
    n = 200
    r = 9.0 # basically duplicates final frame of segment 1
    dtau = 0.0731 # lands us very close to r=3 on final frame
  if segment==3:
    # r=3 to 0.5
    n = 200
    r = 3.0
    dtau = 0.01626
  if segment==4:
    # r=0.5 to 0
    n = 200
    r = 0.5
    dtau = 0.0011918
  for i in range(n):
    outfile = "animation/seg"+("%01d" % segment)+"frame"+("%03d" % i)+".png"
    do_it = (prep_level==2 and (i==0 or i==n-1)) or \
            (prep_level==3 and (i==0 or i==n-1 or i%10==0)) or \
            prep_level==4
    if do_it:
      print("---------------------- r=",r,", file=",outfile," --------------- ",\
           datetime.datetime.now().strftime('%H:%M:%S'))
      sys.stdout.flush()
      # do_image(r,outfile,if_fake,star_catalog_max_mag,width,height,fov_deg,view_rot_deg) # qwe --uncomment
    spacetime = SP_SCH
    chart = CH_SCH
    pars = {}
    x_obs,v_obs,rho,j = ray.schwarzschild_standard_observer(r,spacetime,chart,pars)
    r = r+dtau*v_obs[1]
    if r<0.0:
      break
