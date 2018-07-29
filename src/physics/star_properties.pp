import random

#include "language.h"
#include "math.h"

import math_util

def list_of_spectral_classes():
  return ['O','B','A','F','G','K','M']

def frequency_of_spectral_classes():
  """
  Returns a hash with keys "O", "B", ..., and values that are normalized to 1.
  This is based on the 10,000 brightest stars in the kstars database. This is an estimate of the frequency
  at a given apparent magnitude, and is totally different from the actual frequency of these types in
  the population of stars in the galaxy.
  """
  # The following are minimally rounded, in order to keep from messing up normalization.
  return {"O":0.005417335473515248, "B":0.18148073836276082, "A":0.20676163723916532, "F":0.13764044943820225,\
          "G":0.12981540930979132, "K":0.28230337078651685, "M":0.056581059390048156 }

def random_spectral_class():
  """
  See frequency_of_spectral_class() for how the probabilities are defined.
  """
  freq = frequency_of_spectral_classes()
  names = list_of_spectral_classes()
  n = len(names)
  u = random.random()
  for i in range(n):
    name = names[i]
    f = freq[name]
    if u>f:
      u = u-f
    else:
      return name
  return 'M' # in case normalization fails due to rounding

def spectral_class_to_bv(cl):
  """
  Given a spectral class such as "G", returns the approximate B-V color.
  """
  # https://en.wikipedia.org/wiki/Color_index
  return {'O':-0.33,'B':-0.30,'A':-0.02,'F':0.30,'G':0.58,'K':0.81,'M':1.40}[cl]

def bv_to_log_temperature(bv):
  """
  Given a B-V color, find an approximate blackbody temperature.
  """
  table = [[-0.33,10.65],[-0.30,10.31],[-0.02,9.19],[0.30,8.90],[0.58,8.69],[0.81,8.55],[1.40,8.25]]
  return math_util.linear_interp_from_table(table,0,1,bv,0,len(table)-1)


def log_temperature_to_hue_and_sat(log_temp):
  """
  Given a blackbody temperature, find an approximate hue and saturation. Hue is
  on a scale from 0 (red) to 1 (blue), corresponding to 0-270 in the hsluv system.
  Saturation is from 0 to 1, corresponding to 0-100 in hsluv.
  """
  # subjective tabulation based on:
  #   https://en.wikipedia.org/wiki/File:TernaryColorTmap.PNG
  #   arcturus --  http://www.allthesky.com/constellations/bootes/constell.html
  #   betelgeuse and rigel -- http://www.allthesky.com/constellations/orion/constell.html
  # cf. http://www.tannerhelland.com/4435/convert-temperature-rgb-algorithm-code/
  # The points used are:
  #          500 K - hypothetical, would be saturated red
  #   M -  3,600   - betelgeuse
  #   K -  4,300   - arcturus
  #   A -  7,000   - canopus
  #   B - 12,000   - rigel
  #   O - 38,000   - meissa
  table = [[6.21,0.0,1.0],[8.19,0.13,0.80],[8.37,0.16,0.45],[8.85,0.48,0.30],[9.39,0.81,0.56],[10.55,0.97,0.82]]
  log_temp = math_util.force_into_range(log_temp,6.2,10.8)
  # Force it into a range corresponding to 500 to 50,000 K.
  # At temps below 500 K, it's basically just saturated red.
  # At temps above 50,000 K, you just approach a certain bluish color (it does *not* shift into the violet).
  n = len(table)
  hue = math_util.linear_interp_from_table(table,0,1,log_temp,0,n-1)
  sat = math_util.linear_interp_from_table(table,0,2,log_temp,0,n-1)
  return [hue,sat]

def bv_to_color(bv):
  return log_temperature_to_hue_and_sat(bv_to_log_temperature(bv))
