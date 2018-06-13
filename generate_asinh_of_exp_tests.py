#!/usr/bin/python3

# Generate test values for asinh_of_exp, to be cut and pasted into test_math_util.pp.
# Remove final comma.

# apt-get install mpmath

from mpmath import mp
mp.dps = 50 # decimal precision

def print_test(i):
  x = mp.exp(i)
  y = mp.log(x+mp.sqrt(x**2+1))-i-mp.log(2)
  print("  [",i,",",y,"], # asinh(exp(",i,"))-",i,"-ln(2)=...")

for i in range(1,61,5):
  print_test(i)

