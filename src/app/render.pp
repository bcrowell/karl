#!/usr/bin/python3

#include "language.h"
#include "math.h"

import csv,ephem,json
import euclidean

def main():
  csv_file = "stars.csv"
  outfile = "stars.json"
  verbosity = 1
  #----
  w = 1200 # pixels, width of square image
  h = 600
  fov_deg = 90.0 # horizontal field of view in degrees
  view_rot_deg = 180.0 # 0 means looking at black hole, 180 means looking directly away; positive is pan to right
  #----
  blur = 1 # std dev. of gaussian blur, in units of pixels
  truncate_blur = 4.0 # # truncate blur at 3 s.d.
  exposure = 3.0
  gamma = 0.75 # https://en.wikipedia.org/wiki/Gamma_correction
  # A value less than 1 makes stars appear more uniform in brightness, makes
  # dim stars easier to see.
  #----
  if verbosity>=1:
    PRINT("starting rendering (render.py)")
  #----
  fov_rad = euclidean.deg_to_rad(fov_deg)
  view_rot_rad = euclidean.deg_to_rad(view_rot_deg)
  rot = euclidean.rotation_matrix_from_axis_and_angle(MATH_PI-view_rot_rad,[0.0,1.0,0.0])
  k_proj = (w/2.0)/tan(fov_rad/2.0) # proportionality constant for stereographic projection
  image_i = []
  image_h = []
  image_s = []
  for i in range(w):
    image_i.append(EMPTY1DIM(h)) # intensity
    image_h.append(EMPTY1DIM(h)) # intensity-averaged hue
    image_s.append(EMPTY1DIM(h)) # intensity-averaged saturation
  for i in range(w):
    for j in range(h):
      image_i[i][j] = 0.0
      image_h[i][j] = 0.0
      image_s[i][j] = 0.0
  table = read_csv_file(csv_file)
  for i in range(len(table)):
    row = table[i]
    alpha,phi,raw_brightness,bv = float(row[0]),float(row[1]),float(row[2]),float(row[3])
    if TRUE:
      v = euclidean.spherical_to_cartesian(alpha,phi)
      v = euclidean.apply_matrix(rot,v)
      theta,phi2 = euclidean.cartesian_to_spherical(v)
      r = k_proj*tan(theta/2.0) # stereographic projection, so black hole's silhouette always appears circular
    else:
      r = w*alpha/(MATH_PI/2.0)
      # ... rough and ready projection; away from b.h.; gives roughly 45-degree radius field of view
      phi2 = phi
    x = FLOOR(0.5*w+r*cos(phi2)+0.5)
    y = FLOOR(0.5*h+r*sin(phi2)+0.5)
    if not (x>=0 and x<=w-1 and y>=0 and y<=w-1):
      continue
    brightness = raw_brightness*exposure
    if brightness<=1.0:
      bloat = 1.0
    else:
      bloat = sqrt(brightness) # increase radius of blob to represent greater brightness
    p = CEIL(bloat*blur*truncate_blur+0.5) # size of square and inscribed circle on which to do computations
    blur2 = blur*blur
    for ii in range(2*p+1):
      i = ii-p
      xx = x+i
      if not (xx>=0 and xx<=w-1):
        continue
      for jj in range(2*p+1):
        j = jj-p
        yy = y+j
        if not (yy>=0 and yy<=h-1):
          continue
        # gaussian blur
        r2 = i*i+j*j
        gaussian_x = r2/blur2
        if bloat>1.0:
          gaussian_x=gaussian_x-(bloat-1.0)
          if gaussian_x<0.0:
            gaussian_x = 0.0
        if not gaussian_x<truncate_blur*truncate_blur:
          continue # cut off to a circle, so it doesn't look boxy
        b = brightness*exp(-0.5*gaussian_x)
        image_i[xx][yy] = image_i[xx][yy]+b
        hue,sat = bv_to_color(bv)
        image_h[xx][yy] = image_h[xx][yy]+b*hue
        image_s[xx][yy] = image_s[xx][yy]+b*sat
  max = 0.0
  for i in range(w):
    for j in range(y):
      if image_i[i][j]>max:
        max = image_i[i][j]
  if max==0.0:
    print("max=0")
    exit(-1)
  # Make hue and saturation into intensity-weighted averages:
  for i in range(w):
    for j in range(h):
      z = image_i[i][j]
      if z>0.0:
        image_h[i][j] = image_h[i][j]/z
        image_s[i][j] = image_s[i][j]/z
  # Gamma correction:
  for i in range(w):
    for j in range(h):
      image_i[i][j] = exposure*max*(image_i[i][j]/max)**gamma
  with open(outfile, 'w') as f:
    f.write(json.dumps([image_i,image_h,image_s]))
  if verbosity>=1:
    PRINT("done rendering (render.py)")

def bv_to_color(bv):
  # https://en.wikipedia.org/wiki/Color_index
  # The following is just a rough-and-ready approximation that I made up.
  # Better: http://www.tannerhelland.com/4435/convert-temperature-rgb-algorithm-code/
  # Return values:
  #   hue: 0=red, 1=blue
  #   sat: 0 to 1
  h = (1.40-bv)/(1.40-(-0.33))
  h = put_in_range(h,0,1)
  s = (bv-0.5)**4
  s = put_in_range(s,0,1)
  return [h,s]

def put_in_range(x,min,max):
  if x<min:
    return min
  if x>max:
    return max
  return x

def read_csv_file(filename):
  with open(filename, 'r') as f:
    data = []
    reader = csv.reader(f)
    for row in reader:
      data.append(row)
    return data

main()

