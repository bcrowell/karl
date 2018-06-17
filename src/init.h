# This should only be included (and executed) once, by the main program.

#if "LANG" eq "python"
# Make numpy die with an error on any arithmetic exception, don't just print a warning.
numpy.seterr(all='raise') 
#endif


