      /*
               --- main program test_kruskal ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
        load("lib/loader.js");
      }
      if (typeof test_kruskal === 'undefined') {
        var test_kruskal = {};
      }
      var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps, i, j, k, t0, r0, theta, phi, v, duration, q;

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
      /* The SP_ labels tell us what spacetime we're in. */
      /* The CH_ labels refer to charts within that particular spacetime. */
      /* These are designed so that we can bitwise or them. */
      /* The physics code is written in python, and the js version is automatically translated */
      /* from python, so it has already had these constants substituted in via filepp. But */
      /* For browser-based user interface code written in js, these constants are also */
      /* defined in util/constants.js. */
      /* There is also a spacetimes_c.h version of this file for C sources. */
      /* ... Schwarzschild spacetime */
      /* ... sch5 coordinates */
      /* ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...) */
      /* ... ``Keplerian'' coordinates (t,u,...), with u=r^3/2 */
      /* return codes for Runge-Kutta, designed to be bitwise or-able. */
      /* ... something went really wrong, output is garbage */
      /* ... the geodesic was incomplete */
      /* ... relative precision for arithmetic */
      /* ... ln of (1.0e-16) */
      /* ... ln of greatest number we can store in floating point; IEEE-754 floating point can store 2^128-1 */
      karl.load("kruskal");
      karl.load("transform");
      karl.load("runge_kutta");
      karl.load("angular");
      karl.load("vector");
      karl.load("lambert_w_stuff");
      test_kruskal.test_round_trip_ksk = function(a, b) {
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
        test.assert_rel_equal_eps(a, a2, 10 * (1.0e-16));
        test.assert_rel_equal_eps(b, b2, 10 * (1.0e-16));
      };
      test_kruskal.test_round_trip_sks = function(t, r) {
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
        test.assert_rel_equal_eps(t, t2, 10 * (1.0e-16));
        test.assert_rel_equal_eps(r, r2, 10 * (1.0e-16));
      };
      test_kruskal.simple_free_fall = function() {
        var t0, r0, lambda_max, a, b, x, jac, tdot, v, n, opt, err, final_x, final_v, final_a, final_lambda, info, tf, rf, delta_r, m, g, delta_r_newtonian;

        /* Free fall from rest in a semi-newtonian region, where all the math can be evaluated without */
        /* special massaging to avoid overflows. */
        t0 = 0.0;
        r0 = 300.0; /* as big as possible without causing overflows */
        lambda_max = ((0.01) * (Math.pow((r0), (1.5)))); /* short compared to the Newtonian time required to hit the singularity */
        (function() {
          var temp = transform.schwarzschild_to_kruskal(t0, r0);
          a = temp[0];
          b = temp[1]
        })();
        x = [a, b, 1.0, 0.0, 0.0];
        jac = transform.jacobian_schwarzschild_to_kruskal(t0, r0);
        tdot = 1.0; /* approximate newtonian value of dt/dtau */
        v = [jac[0][0] * tdot, jac[1][0] * tdot, 0.0, 0.0, 0.0];
        n = 100;
        opt = {
          'lambda_max': lambda_max,
          'dlambda': lambda_max / n,
          'ndebug': 0
        };
        (function() {
          var temp = runge_kutta.trajectory_simple(256, 2, x, v, opt);
          err = temp[0];
          final_x = temp[1];
          final_v = temp[2];
          final_a = temp[3];
          final_lambda = temp[4];
          info = temp[5]
        })();
        if (err & 1) {
          throw 'error: ' + info['message'];;
        }
        (function() {
          var temp = transform.kruskal_to_schwarzschild(final_x[0], final_x[1]);
          tf = temp[0];
          rf = temp[1]
        })();
        delta_r = r0 - rf;
        m = 0.5; /* coordinates are such that mass=1/2 */
        g = ((m) * (Math.pow((r0), (-2.0)))); /* approximate as constant accel */
        delta_r_newtonian = ((0.5) * (g) * (Math.pow((lambda_max), (2.0))));
        test.assert_rel_equal_eps(delta_r, delta_r_newtonian, 0.02);
      };
      test_kruskal.test_christoffel_raw_vs_massaged = function(a, b) {
        var p, ch, ch2;

        p = [a, b, 0.8, 0.6, 0.0];
        ch = kruskal.christoffel_raw_maxima_output(p);
        ch2 = kruskal.christoffel_massaged_maxima_output(p);
        for (var i = 0; i < 5; i++) {
          for (var j = 0; j < 5; j++) {
            for (var k = 0; k < 5; k++) {
              /*print("i=",i,", j=",j,", k=",k,", raw=",ch[i][j][k],", massaged=",ch2[i][j][k]) */
              test.assert_rel_equal_eps(ch[i][j][k], ch2[i][j][k], 10 * (1.0e-16));
            }
          }
        }
      };
      test_kruskal.test_motion_kruskal_vs_schwarzschild = function(t0, r0, flip, theta, phi, v, duration) {
        var i, j, k, lambda_max, a0, b0, x0s, x0k, v0s, jac, v0k, n, ndebug, opt, x, v, chart, err, final_x, final_v, final_a, final_lambda, info, xfs, xfk, tf, rf, eps;

        /* initial point = t0,r0,theta,phi */
        /* initial velocity = v, expressed in (t,r,i,j,k) component form */
        /* If flip is true, then the Kruskal version goes to (-a,-b). */
        /* For convenience, v will automatically be normalized and made tangent to the 2-sphere. */
        /* The variable duration is the maximum affine parameter, in units of the newtonian */
        /* free-fall time required in order to hit the singularity from this radius. */
        (function() {
          var temp = angular.theta_phi_to_ijk(theta, phi);
          i = temp[0];
          j = temp[1];
          k = temp[2]
        })();
        lambda_max = ((duration) * (Math.pow((r0), (1.5))));
        (function() {
          var temp = transform.schwarzschild_to_kruskal(t0, r0);
          a0 = temp[0];
          b0 = temp[1]
        })();
        if (flip) {
          a0 = -a0;
          b0 = -b0;
        }
        /* --- Find initial position in both coordinate systems: ------------------------------------ */
        x0s = (karl.clone_array1d(([t0, r0, i, j, k]))); /* initial point in Schwarzschild coordinates */
        x0k = (karl.clone_array1d(([a0, b0, i, j, k]))); /* ... in Kruskal */
        /* --- Find initial velocity vector in Schwarzschild coordinates: --------------------------- */
        v0s = angular.make_tangent(x0s, v); /* initial velocity in Schwarzschild coordinates */
        v0s = vector.normalize(256, 1, x0s, v0s);
        /* --- Find initial velocity vector in Kruskal coordinates: --------------------------------- */
        jac = transform.jacobian_schwarzschild_to_kruskal(t0, r0);
        v0k = (karl.clone_array1d(v0s));
        v0k[0] = jac[0][0] * v0s[0] + jac[0][1] * v0s[1];
        v0k[1] = jac[1][0] * v0s[0] + jac[1][1] * v0s[1];
        v0k = vector.normalize(256, 2, x0k, v0k);
        if (flip) {
          v0k[0] = -v0k[0];
          v0k[1] = -v0k[1];
        }
        /* --- */
        if (verbosity >= 3) {
          print("----------------------");
          print("x0s=", io_util.vector_to_str_n_decimals(x0s, 8));
          print("x0k=", io_util.vector_to_str_n_decimals(x0k, 8));
          print("v0s=", io_util.vector_to_str_n_decimals(v0s, 8));
          print("v0k=", io_util.vector_to_str_n_decimals(v0k, 8));
        }
        /* --- */
        for (var i = 0; i < 2; i++) {
          n = 100; /* fails with n=1000, why? */
          ndebug = 0;
          if ((false) && verbosity >= 3) {
            print("----------------------");
            ndebug = n / 100;
          }
          opt = {
            'lambda_max': lambda_max,
            'dlambda': lambda_max / n,
            'ndebug': ndebug
          };
          if (i == 0) {
            x = x0s;
            v = v0s;
            chart = 1;
          } else {
            x = x0k;
            v = v0k;
            chart = 2;
          }
          (function() {
            var temp = runge_kutta.trajectory_simple(256, chart, x, v, opt);
            err = temp[0];
            final_x = temp[1];
            final_v = temp[2];
            final_a = temp[3];
            final_lambda = temp[4];
            info = temp[5]
          })();
          if (err & 1) {
            throw 'error: ' + info['message'];;
          }
          if (i == 0) {
            xfs = (karl.clone_array1d(final_x));
          } else {
            xfk = (karl.clone_array1d(final_x));
          }
        }
        /* Convert final result of Kruskal calculations to Schwarzschild coordinates: */
        (function() {
          var temp = transform.kruskal_to_schwarzschild(xfk[0], xfk[1]);
          tf = temp[0];
          rf = temp[1]
        })();
        xfk[0] = tf;
        xfk[1] = rf;
        if (verbosity >= 3) {
          print("xfs=", io_util.vector_to_str_n_decimals(xfs, 8));
          print("xfk=", io_util.vector_to_str_n_decimals(xfk, 8));
        }
        eps = ((1000.0) * (Math.pow((n), (-4.0))));
        test.assert_rel_equal_eps_vector(xfs, xfk, eps);
      };
      /*-------------------------------------------------------------------------- */
      /* run tests */
      /*-------------------------------------------------------------------------- */
      test_kruskal.test_christoffel_symmetry = function(ch) {
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
      (function() {
        var temp = angular.theta_phi_to_ijk(0.333, 0.1776);
        i = temp[0];
        j = temp[1];
        k = temp[2]
      })();
      test_kruskal.test_christoffel_symmetry(kruskal.christoffel([1.0, -3.3, i, j, k])); /* check symmetry at a random point */
      /* region II: */
      test_kruskal.test_round_trip_ksk(0.5, 0.5);
      /* ... Note that if we test this with coordinates like (a,b)=(epsilon,epsilon), we get */
      /*     poor relative precision, because in Schwarzschild coordinates, r=1+epsilon. */
      /* region I: */
      test_kruskal.test_round_trip_sks(0.0, 2.0);
      test_kruskal.test_round_trip_ksk(1.0, -1.0);
      test_kruskal.test_round_trip_sks(100.0, 100.0);
      test_kruskal.test_christoffel_raw_vs_massaged(0.5, -0.5); /* region I */
      test_kruskal.test_christoffel_raw_vs_massaged(0.5, 0.5); /* region II */
      test_kruskal.test_christoffel_raw_vs_massaged(-0.5, 0.5); /* region III */
      test_kruskal.test_christoffel_raw_vs_massaged(-0.5, -0.5); /* region IV */
      test_kruskal.simple_free_fall();
      /*--------------------------------------------------------------- */
      /* a bunch of applications of test_motion_kruskal_vs_schwarzschild: */
      /*--------------------------------------------------------------- */
      t0 = 0.0;
      r0 = 2.0;
      theta = Math.PI / 2.0;
      phi = 0.0;
      v = [1.0, 0.0, 0.0, 0.0, 0.0]; /* initially at rest; this gets normalized later */
      duration = 0.2;
      test_kruskal.test_motion_kruskal_vs_schwarzschild(t0, r0, (false), theta, phi, v, duration); /* region I */
      test_kruskal.test_motion_kruskal_vs_schwarzschild(t0, r0, (true), theta, phi, v, duration); /* region III */
      r0 = 0.5;
      v = [0.0, -1.0, 0.0, 0.0, 0.0];
      test_kruskal.test_motion_kruskal_vs_schwarzschild(t0, r0, (false), theta, phi, v, duration); /* region II */
      test_kruskal.test_motion_kruskal_vs_schwarzschild(t0, r0, (true), theta, phi, v, duration); /* region IV */
      /* A test with random values of everything, to exercise all Christoffel symbols: */
      r0 = 2.0;
      theta = 1.111;
      phi = 2.345;
      v = [1.0, 0.01776, 0.01066, 0.01492, 0.02001];
      duration = 0.2;
      test_kruskal.test_motion_kruskal_vs_schwarzschild(t0, r0, (false), theta, phi, v, duration);
      /* Large r and t, test whether we get any overflows: */
      t0 = 0.0;
      r0 = 1.0e8;
      theta = Math.PI / 2.0;
      phi = 0.0;
      q = 1 / r0; /* scale angular motion down by this amount to keep the motion from being FTL */
      /*v = [1.0,0.01776,0.01066*q,0.01492*q,0.02001*q] # random initial motion */
      v = [1.0, 0, 0, 0, 0]; /* random initial motion */
      duration = 1.0e-10; /* Kruskal coordinates are not well adapted to covering the motion of */
      /*                    a nonrelativistic object for long times. */
      test_kruskal.test_motion_kruskal_vs_schwarzschild(t0, r0, (false), theta, phi, v, duration);
      test.done(verbosity, "kruskal");
