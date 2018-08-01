argv = ["js/test_angular.jsi", "karl"]
/*
   --- main program test_angular ---
   This was translated from python. Do not edit directly.
*/
var karl = {};
karl.modules_loaded = {};
if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
  load("lib/loader.js");
}
if (typeof test_angular === 'undefined') {
  var test_angular = {};
}
var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps, x, n, v;

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
karl.load("angular");
x = [0, 0, 1, 2, 3];
n = angular.renormalize(x);
test.assert_equal((((Math.pow((n[2]), (2.0))) + (Math.pow((n[3]), (2.0))) + (Math.pow((n[4]), (2.0))))), (1.0));
x = [0, 0, 1, 0, 0];
v = [0, 0, 3, 700, 0];
v = angular.make_tangent(x, v);
test.assert_equal(v[2], 0);
test.assert_equal(v[3], 700);
test.done(verbosity, "angular");
