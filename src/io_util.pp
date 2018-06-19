import sys

# l is an array of objects which may not be strings
def strcat(l):
  return ''.join(map(str, l)) \
  #js return l.join("");

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

#if "LANG" eq "python"
def vector_to_str_n_decimals(v,n):
  f = []
  for x in v:
    fmt = "%"+str(n+2)+"."+str(n)+"e"
    f.append(fmt % x)
  return ','.join(f)
#endif
#js function vector_to_str_n_decimals(v,n) {
#js   var s = "";
#js   for (var i=0; i<v.length; i++) {
#js     s = s+v[i].toExponential(n).toString();
#js   }
#js   return s;
#js }

def fl(x):
  return "%5.3e" % x #js return x.toExponential(3).toString();
