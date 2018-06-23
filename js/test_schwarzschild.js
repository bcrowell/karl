      /*
               --- main program test_schwarzschild ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
          load("lib/loader.js");
      }
      if (typeof test_schwarzschild === 'undefined') {
          var test_schwarzschild = {};
      }
      var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps, t, r, sin_theta, cos_theta, i, j, k, ch4, ch5;

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
      karl.load("schwarzschild");
      test_schwarzschild.test_christoffel_symmetry = function(ch) {
          var l;

          l = ch.length;
          for (var i = 0; i < l; i++) {
              for (var j = 0; j < l; j++) {
                  for (var k = 0; k < l; k++) {
                      test.assert_equal(ch[i][j][k], ch[j][i][k]);
                  }
              }
          }
      };
      /* Exercise as many lines of code as possible, using random coordinates. */
      t = 3.3;
      r = 2.2;
      sin_theta = 1 / Math.sqrt(2.0);
      cos_theta = 1 / Math.sqrt(2.0);
      /* A set of i,j,k that lie on the unit sphere: */
      i = 0.6;
      j = 0.8;
      k = 0.0;
      schwarzschild.metric_sch4(r, sin_theta);
      schwarzschild.sigma(r);
      ch4 = schwarzschild.christoffel_sch4(t, r, sin_theta, cos_theta);
      ch5 = schwarzschild.christoffel([t, r, i, j, k]);
      /* test symmetry of Christoffel symbols */
      test_schwarzschild.test_christoffel_symmetry(ch4);
      test_schwarzschild.test_christoffel_symmetry(ch5);
      /* As far as testing whether the results are actually correct, see */
      /* tests in test_runge_kutta. */
      test.done(verbosity, "schwarzschild");
