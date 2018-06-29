// The SP_ labels tell us what spacetime we're in.
// The CH_ labels refer to charts within that particular spacetime.
// These are designed so that we can bitwise or them.
// The physics code is written in python, and the js version is automatically translated
// from python, so it has already had these constants substituted in via filepp, which
// gets them from spacetimes.h

var SP_SCH=256; // Schwarzschild spacetime
var CH_SCH=1; // sch5 coordinates
var CH_AKS=2; // Kruskal-Szekeres null coordinates (asinh V,asinh W,...),

