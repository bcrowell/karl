/*
   Code snippets extracted from
     https://github.com/josdejong/mathjs ,
   Jos de Jong.
   Apache License, Version 2.0
*/
if (typeof(Math.cosh) === 'undefined') {
  /* ECMAScript 6 defines hyperbolic trig functions. */
  Math.cosh = function(x) { return (Math.exp(x)+Math.exp(-x))/2};
  Math.sinh = function(x) { return (Math.exp(x)-Math.exp(-x))/2};
  Math.tanh = function(x) { return (Math.exp(x)-Math.exp(-x))/(Math.exp(x)+Math.exp(-x))};
}
if (typeof(Math.acosh) === 'undefined') {
  /* ECMAScript 6 defines inverse hyperbolic trig functions. */
  Math.arccosh = function(x) { return Math.log(Math.sqrt(x * x - 1) + x) };
  Math.arcsinh = function(x) { return Math.log(Math.sqrt(x * x + 1) + x) };
  Math.arctanh = function(x) { return Math.log((1 + x) / (1 - x)) / 2 };
}
/* Add aliases like "arccosh" rather than "acosh" to be consistent with numpy. */
Math.arccosh = Math.acosh;
Math.arcsinh = Math.asinh;
Math.arctanh = Math.atanh;


