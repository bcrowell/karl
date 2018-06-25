# Stuff for calling C subroutines.

#if "LANG" eq "python"

import ctypes,os,inspect

mydir = os.path.dirname(os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename))
# ... the directory in which the currently running script is -- https://stackoverflow.com/a/6209894/1142217
karl_c_lib = ctypes.cdll.LoadLibrary(mydir+'/karl.so')

c_double_p = ctypes.POINTER(ctypes.c_double)

#endif


