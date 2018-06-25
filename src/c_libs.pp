# Stuff for calling C subroutines.

#if "LANG" eq "python"
import ctypes
karl_c_lib = ctypes.cdll.LoadLibrary('/home/bcrowell/Documents/programming/karl/obj/karl.so')
# ... fixme: don't hardcode the directory
c_double_p = ctypes.POINTER(ctypes.c_double)
#endif


