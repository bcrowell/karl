#!/usr/bin/python3

#include "language.h"
#include "math.h"

import csv,ephem,json

def main():
  csv_file = "stars.csv"
  outfile = "stars.json"
  #----
  w = 1000 # pixels, width of square image
  fov_deg = 45.0 # field of view in degrees
  view_rot_deg = 0.0 # angle away from black hole, i.e., 180 means looking directly away from it
  #----
  blur = 1 # std dev. of gaussian blur, in units of pixels
  truncate_blur = 4.0 # # truncate blur at 3 s.d.
  overexpose = 10.0
  gamma = 0.75 # https://en.wikipedia.org/wiki/Gamma_correction
  # A value less than 1 makes stars appear more uniform in brightness, makes
  # dim stars easier to see.
  #----
  image_i = EMPTY2DIM(w) # intensity
  image_h = EMPTY2DIM(w) # intensity-averaged hue
  image_s = EMPTY2DIM(w) # intensity-averaged saturation
  for i in range(w):
    for j in range(w):
      image_i[i][j] = 0.0
      image_h[i][j] = 0.0
      image_s[i][j] = 0.0
  table = read_csv_file(csv_file)
  for i in range(len(table)):
    row = table[i]
    alpha,phi,brightness,bv = float(row[0]),float(row[1]),float(row[2]),float(row[3])
    if FALSE:
      r = w*alpha/(MATH_PI/2.0)
      # ... rough and ready projection; away from b.h.; gives roughly 45-degree radius field of view
    else:
      r = w*(MATH_PI-alpha)/(MATH_PI)
      # towards b.h., bigger field of view
    x = FLOOR(0.5*w+r*cos(phi)+0.5)
    y = FLOOR(0.5*w-r*sin(phi)+0.5)
    if not (x>=0 and x<=w-1 and y>=0 and y<=w-1):
      continue
    p = CEIL(blur*truncate_blur+0.5) # size of square and inscribed circle on which to do computations
    blur2 = blur*blur
    for ii in range(2*p+1):
      i = ii-p
      xx = x+i
      if not (xx>=0 and xx<=w-1):
        continue
      for jj in range(2*p+1):
        j = jj-p
        yy = y+j
        if not (yy>=0 and yy<=w-1):
          continue
        # gaussian blur
        h = i*i+j*j
        if not h/blur2<truncate_blur*truncate_blur:
          continue # cut off to a circle, so it doesn't look boxy
        b = brightness*exp(-0.5*h/blur2)
        image_i[xx][yy] = image_i[xx][yy]+b
        hue,sat = bv_to_color(bv)
        image_h[xx][yy] = image_h[xx][yy]+b*hue
        image_s[xx][yy] = image_s[xx][yy]+b*sat
  max = 0.0
  for i in range(w):
    for j in range(w):
      if image_i[i][j]>max:
        max = image_i[i][j]
  if max==0.0:
    print("max=0")
    exit(-1)
  # Make hue and saturation into intensity-weighted averages:
  for i in range(w):
    if math.isinf(image_h[487][504]):
      THROW('qwe250')
    for j in range(w):
      z = image_i[i][j]
      if z>0.0:
        image_h[i][j] = image_h[i][j]/z
        image_s[i][j] = image_s[i][j]/z
  # Gamma correction:
  for i in range(w):
    for j in range(w):
      image_i[i][j] = max*(image_i[i][j]/max)**gamma
  with open(outfile, 'w') as f:
    f.write(json.dumps([image_i,image_h,image_s]))

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

