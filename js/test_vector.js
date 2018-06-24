      /*
               --- main program test_vector ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
          load("lib/loader.js");
      }
      if (typeof test_vector === 'undefined') {
          var test_vector = {};
      }
      var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps, d, max_err;

      /*!/usr/bin/python3 */
      /* ... note that (NaN)==(NaN) is false, so use IS_(NaN) */
      karl.load("lib/array");;
      /* ... works in rhino and d8 */
      /* ... https://stackoverflow.com/q/26738943/1142217 */
      /*           ... see notes above about usage with array literals */
      /*           ... in JS, numbers are primitives, not objects, so no need clone them */
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
      /* The SP_ labels tell us what spacetime we're in. */
      /* The CH_ labels refer to charts within that particular spacetime. */
      /* These are designed so that we can bitwise or them. */
      /* The physics code is written in python, and the js version is automatically translated */
      /* from python, so it has already had these constants substituted in via filepp. But */
      /* For browser-based user interface code written in js, these constants are also */
      /* defined in util/constants.js. */
      /* ... Schwarzschild spacetime */
      /* ... sch5 coordinates */
      /* ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...),  */
      /* ... relative precision for arithmetic */
      /* ... ln of (1.0e-16) */
      /* ... ln of greatest number we can store in floating point; IEEE-754 floating point can store 2^128-1 */
      karl.load("vector");
      karl.load("transform");
      test_vector.test_norm_schwarzschild_vs_kruskal = function(t, r, dt, dr, eps) {
          var p_s, norm_s, a, b, a2, b2, p_k, norm_k;

          p_s = [t, r, 0, 0]; /* point in Schwarzschild coords */
          norm_s = vector.norm4(256, 1, p_s, [dt, dr, 0, 0]);
          /* Transform to arcsinh-Kruskal and compare norms: */
          (function() {
              var temp = transform.schwarzschild_to_kruskal(t, r);
              a = temp[0];
              b = temp[1]
          })();
          (function() {
              var temp = transform.schwarzschild_to_kruskal(t + dt, r + dr);
              a2 = temp[0];
              b2 = temp[1]
          })();
          p_k = [a, b, 0, 0];
          norm_k = vector.norm4(256, 2, p_k, [a2 - a, b2 - b, 0, 0]);
          if (verbosity >= 3) {
              print("test_norm_schwarzschild_vs_kruskal: t=", t, ", r=", r, ", a=", a, ", b=", b, ", norm_s=", norm_s, ", norm_k=", norm_k);
          }
          test.assert_rel_equal_eps(norm_s, norm_k, eps);
      };
      d = Math.sqrt((1.0e-16));
      max_err = 20 * d;
      test_vector.test_norm_schwarzschild_vs_kruskal(0.111, 2.0, d, 0.777 * d, max_err); /* random point and displacement in region I */
      test_vector.test_norm_schwarzschild_vs_kruskal(0.111, 0.5, d, 0.777 * d, max_err); /* ... region II */
      test.done(verbosity, "vector");
