# The SP_ labels tell us what spacetime we're in.
# The CH_ labels refer to charts within that particular spacetime.
# These are designed so that we can bitwise or them.
# The physics code is written in python, and the js version is automatically translated
# from python, so it has already had these constants substituted in via filepp. But
# For browser-based user interface code written in js, these constants are also
# defined in util/constants.js.
# There is also a spacetimes_c.h version of this file for C sources.

#define SP_SCH 256
# ... Schwarzschild spacetime

#define CH_SCH 1
# ... sch5 coordinates
#define CH_AKS 2
# ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...), 
