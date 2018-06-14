# The SP_ labels tell us what spacetime we're in.
# The CH_ labels refer to charts within that particular spacetime.
# These are designed so that we can bitwise or them.

#define SP_SCH 256

#define CH_SCH 1
  # ... sch5 coordinates
#define CH_AKS 2
  # ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...), 
