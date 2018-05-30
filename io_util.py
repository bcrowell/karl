import sys

# l is an array of objects which may not be strings
def strcat(l):
  return ''.join(map(str, l))

def print_no_newline(s):
  sys.stdout.write(s)
  sys.stdout.flush()

def vector_to_str(v):
  f = []
  for x in v:
    f.append("%5.3e" % x)
  return ','.join(f)
