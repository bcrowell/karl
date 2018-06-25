# Stuff for calling C subroutines.

#if "LANG" eq "python"

import ctypes,os,inspect

mydir = os.path.dirname(os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename))
# ... the directory in which the currently running script is -- https://stackoverflow.com/a/6209894/1142217
karl_c_lib = ctypes.cdll.LoadLibrary(mydir+'/karl.so')

c_double_p = ctypes.POINTER(ctypes.c_double)

karl_c_lib.veberic_lambert_w.restype = ctypes.c_double
karl_c_lib.veberic_lambert_w.argtypes = [ctypes.c_double]
karl_c_lib.lambert_w_of_exp.restype = ctypes.c_double
karl_c_lib.lambert_w_of_exp.argtypes = [ctypes.c_double]
karl_c_lib.aux.argtypes = [c_double_p,c_double_p,c_double_p,ctypes.c_double,ctypes.c_double]
#endif


