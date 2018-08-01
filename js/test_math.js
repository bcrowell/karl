argv = ["js/test_math.jsi", "karl"]
/*
   --- main program test_math ---
   This was translated from python. Do not edit directly.
*/
var karl = {};
karl.modules_loaded = {};
if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
  load("lib/loader.js");
}
if (typeof test_math === 'undefined') {
  var test_math = {};
}
var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps;

/*!/usr/bin/python3 */
/* ... note that (NaN)==(NaN) is false, so use IS_(NaN) */
karl.load("lib/array");;
/* ... works in rhino and d8 */
/* ... https://stackoverflow.com/q/26738943/1142217 */
/* ... usage: throw io_util.strcat(([...])); ... extra parens required by filepp so it believes it's a single argument */
/*           ... see notes above about usage with array literals */
/*                 ... works in rhino */
if (!(typeof window !== 'undefined') && (typeof Math.karl === 'undefined')) {
  /* load() works in rhino,  !  sure about other engines */
  load("lib/math.js");
  if (typeof one_over_E === 'undefined') {
    load("lib/lambertw.js")
  }
  load("lambert_w_stuff.js");
  Math.karl = 1;
}
/* This should only be included (and executed) once, by the main program. */
load("test.js");;
assert_rel_equal = test.assert_rel_equal;;
assert_equal = test.assert_equal;;
assert_rel_equal_eps = test.assert_rel_equal_eps;;
assert_equal_eps = test.assert_equal_eps;;
/* ... relative precision for arithmetic */
/* ... ln of (1.0e-16) */
/* ... ln of greatest number we can store in floating point; IEEE-754 floating point can store 2^128-1 */

test.assert_rel_equal_eps(1776, Math.cosh(Math.arccosh(1776)), 5 * (1.0e-16));
test.assert_rel_equal_eps(1776, Math.sinh(Math.arcsinh(1776)), 5 * (1.0e-16));
test.assert_rel_equal_eps(0.1776, Math.tanh(Math.arctanh(0.1776)), 5 * (1.0e-16));
test.done(verbosity, "test_math");
