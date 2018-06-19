/*
   Code snippets extracted from
     https://github.com/josdejong/mathjs ,
   Jos de Jong.
   Apache License, Version 2.0
*/
/* Names like "arccosh" rather than "acosh" are chosen to be consistent with numpy. */
Math.arccosh = function(x) { return Math.log(Math.sqrt(x * x - 1) + x) }
Math.arcsinh = function(x) { return Math.log(Math.sqrt(x * x + 1) + x) }
Math.arctanh = function(x) { return Math.log((1 + x) / (1 - x)) / 2 }

