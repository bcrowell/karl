      /*
               --- main program test_lambert_w ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
        load("lib/loader.js");
      }
      if (typeof test_lambert_w === 'undefined') {
        var test_lambert_w = {};
      }
      var assert_boolean, assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps, z, n_implementations, u, w;

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
      /* ... relative precision for arithmetic */
      /* ... ln of (1.0e-16) */
      /* ... ln of greatest number we can store in floating point; IEEE-754 floating point can store 2^128-1 */
      /* This should only be included (and executed) once, by the main program. */
      load("test.js");;
      assert_boolean = test.assert_boolean;;
      assert_rel_equal = test.assert_rel_equal;;
      assert_equal = test.assert_equal;;
      assert_rel_equal_eps = test.assert_rel_equal_eps;;
      assert_equal_eps = test.assert_equal_eps;;
      karl.load("lambert_w_stuff");
      z = 0.1492;
      test.assert_rel_equal_eps(z, Math.lambert_w(z * Math.exp(z)), 5 * (1.0e-16));
      n_implementations = 1;
      for (var implementation = 0; implementation < n_implementations; implementation++) {
        u = 1.0;
        for (var i = 0; i < 90; i++) {
          if (implementation == 0) {
            w = lambert_w_stuff.lambert_w_of_exp(u);
          }
          if (verbosity >= 3) {
            print("====== In test_lambert_w, implementation=", implementation, " u=", u, ", w=", w);
          }
          assert_rel_equal(u, Math.log(w) + w);
          u = u * 2.5; /* a little less than e, to try to work out all possible whole-number parts of ln(u) */
        }
      };
      test.done(verbosity, "lambert_w");
