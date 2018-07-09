      /*
               --- main program test_transform ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
        load("lib/loader.js");
      }
      if (typeof test_transform === 'undefined') {
        var test_transform = {};
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
      karl.load("kruskal");
      karl.load("transform");
      test_transform.test_round_trip_ksk = function(a, b) {
        var t, r, mu, a2, b2;

        /* We don't expect this to work if we're in regions III or IV, which can't be represented in S coords. */
        (function() {
          var temp = kruskal.aux(a, b);
          t = temp[0];
          r = temp[1];
          mu = temp[2]
        })();
        (function() {
          var temp = transform.schwarzschild_to_kruskal(t, r);
          a2 = temp[0];
          b2 = temp[1]
        })();
        if (verbosity >= 3) {
          print("ksk: a=", a, ", b=", b, ", t=", t, ", r=", r, ", a'=", a2, ", b'=", b2);
        }
        test.assert_rel_equal_eps(a, a2, 5 * (1.0e-16));
        test.assert_rel_equal_eps(b, b2, 5 * (1.0e-16));
      };
      test_transform.test_round_trip_sks = function(t, r) {
        var a, b, t2, r2, mu;

        (function() {
          var temp = transform.schwarzschild_to_kruskal(t, r);
          a = temp[0];
          b = temp[1]
        })();
        (function() {
          var temp = kruskal.aux(a, b);
          t2 = temp[0];
          r2 = temp[1];
          mu = temp[2]
        })();
        if (verbosity >= 3) {
          print("sks: t=", t, ", r=", r, ", a=", a, ", b=", b, ", t2=", t2, ", r2=", r2);
        }
        test.assert_rel_equal_eps(t, t2, 5 * (1.0e-16));
        test.assert_rel_equal_eps(r, r2, 5 * (1.0e-16));
      };
      test_transform.test_jacobian_schwarzschild_to_kruskal = function(t, r) {
        var jac, d, a0, b0, t1, r1, a1, b1, dk, der;

        /* Test analytic Jacobian against numerical differentiation. */
        jac = transform.jacobian_schwarzschild_to_kruskal(t, r);
        d = Math.sqrt((1.0e-16)); /* differences of 10^-8, so that we can get numerical approximations to derivatives to 8 decimals */
        (function() {
          var temp = transform.schwarzschild_to_kruskal(t, r);
          a0 = temp[0];
          b0 = temp[1]
        })();
        for (var i = 0; i < 2; i++) {
          for (var j = 0; j < 2; j++) {
            t1 = t;
            r1 = r;
            if (j == 0) {
              t1 += d;
            }
            if (j == 1) {
              r1 += d;
            }
            (function() {
              var temp = transform.schwarzschild_to_kruskal(t1, r1);
              a1 = temp[0];
              b1 = temp[1]
            })();
            if (i == 0) {
              dk = a1 - a0;
            }
            if (i == 1) {
              dk = b1 - b0;
            }
            der = dk / d; /* numerical approximation to the partial derivative */
            if (verbosity >= 3) {
              print("i=", i, ", j=", j, ", der=", der, ", J=", jac[i][j]);
            }
            test.assert_rel_equal_eps(der, jac[i][j], 10 * d);
          }
        }
      };
      test_transform.test_jacobian_kruskal_to_schwarzschild = function(t, r) {
        var jac, d, a, b, a1, b1, t1, r1, ds, der;

        /* Test analytic Jacobian against numerical differentiation. */
        jac = transform.jacobian_kruskal_to_schwarzschild(t, r);
        d = Math.sqrt((1.0e-16)); /* differences of 10^-8, so that we can get numerical approximations to derivatives to 8 decimals */
        (function() {
          var temp = transform.schwarzschild_to_kruskal(t, r);
          a = temp[0];
          b = temp[1]
        })();
        for (var i = 0; i < 2; i++) {
          for (var j = 0; j < 2; j++) {
            a1 = a;
            b1 = b;
            if (j == 0) {
              a1 += d;
            }
            if (j == 1) {
              b1 += d;
            }
            (function() {
              var temp = transform.kruskal_to_schwarzschild(a1, b1);
              t1 = temp[0];
              r1 = temp[1]
            })();
            if (i == 0) {
              ds = t1 - t;
            }
            if (i == 1) {
              ds = r1 - r;
            }
            der = ds / d; /* numerical approximation to the partial derivative */
            if (verbosity >= 3) {
              print("i=", i, ", j=", j, ", der=", der, ", J=", jac[i][j]);
            }
            test.assert_rel_equal_eps(der, jac[i][j], 10 * d);
          }
        }
      };
      /* region II: */
      test_transform.test_round_trip_ksk(0.5, 0.5);
      /* ... Note that if we test this with coordinates like (a,b)=(epsilon,epsilon), we get */
      /*     poor relative precision, because in Schwarzschild coordinates, r=1+epsilon. */
      /* region I: */
      test_transform.test_round_trip_sks(0.0, 2.0);
      test_transform.test_round_trip_ksk(1.0, -1.0);
      test_transform.test_round_trip_sks(100.0, 100.0);
      /* In the following, the arguments are (t,r). */
      test_transform.test_jacobian_schwarzschild_to_kruskal(0.111, 2.0); /* region I */
      test_transform.test_jacobian_schwarzschild_to_kruskal(0.222, 0.5); /* region II */
      test_transform.test_jacobian_kruskal_to_schwarzschild(0.111, 2.0); /* region I */
      test_transform.test_jacobian_kruskal_to_schwarzschild(0.222, 0.5); /* region II */
      test.done(verbosity, "transform");
