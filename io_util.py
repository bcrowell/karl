import sys

# l is an array of objects which may not be strings
def strcat(l):
  return ''.join(map(str, l))

def print_no_newline(s):
  sys.stdout.write(s)
  sys.stdout.flush()

def vector_to_str(v):
  return vector_to_str_n_decimals(v,3)

def vector_to_str_n_decimals(v,n):
  f = []
  for x in v:
    fmt = "%"+str(n+2)+"."+str(n)+"e"
    f.append(fmt % x)
  return ','.join(f)
