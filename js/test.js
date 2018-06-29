      /*
               --- module test ---
               This was translated from python. Do not edit directly.
            */
      if (typeof test === 'undefined') {
        var test = {};
      }
      var verbosity;

      karl.load("io_util");

      load("io_util.js");
      /* ... note that (NaN)==(NaN) is false, so use IS_(NaN) */
      karl.load("lib/array");;
      /* ... works in rhino and d8 */
      /* ... https://stackoverflow.com/q/26738943/1142217 */
      /* ... usage: throw io_util.strcat(([...])); ... extra parens required by filepp so it believes it's a single argument */
      /*           ... see notes above about usage with array literals */
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
      test.assert_equal_eps = function(x, y, eps) {
        var err;

        err = x - y;
        if ((isNaN(x)) || (isNaN(y)) || Math.abs(err) > eps) {
          throw io_util.strcat((["test failed, x=", x, ", y=", y, ", err=", err, ", eps=", eps]));;
        }
      };
      test.assert_rel_equal_eps = function(x, y, eps) {
        var rel_err;

        if (x == 0.0 && y == 0.0) {
          return;
        }
        if (x == 0.0) {
          return test.assert_rel_equal_eps(y, x, eps); /* avoid division by zero */
        }
        rel_err = (x - y) / x;
        if ((isNaN(x)) || (isNaN(y)) || Math.abs(rel_err) > eps) {
          throw io_util.strcat((["test failed, x=", x, ", y=", y, ", rel err=", rel_err, ", eps=", eps]));;
        }
      };
      test.assert_equal = function(x, y) {

        return test.assert_equal_eps(x, y, 2.0 * (1.0e-16));
      };
      test.assert_rel_equal = function(x, y) {

        return test.assert_rel_equal_eps(x, y, 2.0 * (1.0e-16));
      };
      test.assert_rel_equal_eps_vector = function(x, y, eps) {

        for (var i = 0; i < x.length; i++) {
          test.assert_rel_equal_eps(x[i], y[i], eps);
        }
      };
      test.done = function(verbosity, name) {

        if (verbosity >= 1) {
          print("Passed test_" + name + ", language=js");
        }
      };
      verbosity = 1;
