import sys

# l is an array of objects which may not be strings
def strcat(l):
  return ''.join(map(str, l))

def print_no_newline(s):
#if "LANG" eq "python"
  sys.stdout.write(s)
  sys.stdout.flush()
#endif
#if "LANG" eq "js"
  PRINT(s)
#endif

def vector_to_str(v):
  return vector_to_str_n_decimals(v,3)

def vector_to_str_n_decimals(v,n):
  f = []
  for x in v:
    fmt = "%"+str(n+2)+"."+str(n)+"e"
    f.append(fmt % x)
  return ','.join(f)

def fl(x):
  return "%5.3e" % x
