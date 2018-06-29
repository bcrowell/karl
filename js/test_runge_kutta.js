      /*
               --- main program test_runge_kutta ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
        load("lib/loader.js");
      }
      if (typeof test_runge_kutta === 'undefined') {
        var test_runge_kutta = {};
      }
      var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps;

      /*!/usr/bin/python3 */
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
      /* ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...),  */
      /* return codes for Runge-Kutta, designed to be bitwise or-able. */
      /* ... something went really wrong, output is garbage */
      /* ... the geodesic was incomplete */
      karl.load("runge_kutta");
      karl.load("angular");
      karl.load("vector");
      /*verbosity=2 */
      /*-------------------------------------------------------------------------------------------------- */
      test_runge_kutta.smoke_test = function() {
        var x, v, ndebug, opt, err, final_x, final_v, final_a, final_lambda, info;

        /* Free fall from rest, from 3 Schwarzschild radii.  */
        /* Affine param is approximately but not exactly equal to proper time. */
        x = [0.0, 10.0, 1.0, 0.0, 0.0];
        v = [1.0, 0.0, 0.0, 0.0, 0.0];
        ndebug = 0;
        if (verbosity >= 3) {
          ndebug = 1;
        }
        opt = {
          'lambda_max': 100.0,
          'dlambda': 10.0,
          'ndebug': ndebug
        };
        (function() {
          var temp = runge_kutta.geodesic_simple(256, 1, x, v, opt);
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
        /* Angular coordinates shouldn't have changed: */
        test.assert_equal(x[2], final_x[2]);
        test.assert_equal(x[3], final_x[3]);
        test.assert_equal(x[4], final_x[4]);
      };
      /*-------------------------------------------------------------------------------------------------- */
      test_runge_kutta.simple_newtonian_free_fall = function() {
        var r0, lambda_max, x, v, n, opt, err, final_x, final_v, final_a, final_lambda, info, rf, delta_r, m, g, delta_r_newtonian, rel_err, expect_err, v_newtonian, vf, rel_err_v;

        /* Newtonian limit. Free fall from ~1 a.u. */
        r0 = 1.0e8;
        lambda_max = ((0.01) * (Math.pow((r0), (1.5)))); /* short compared to the time required to hit the singularity */
        x = [0.0, r0, 1.0, 0.0, 0.0];
        v = [1.0, 0.0, 0.0, 0.0, 0.0];
        n = 100;
        opt = {
          'lambda_max': lambda_max,
          'dlambda': lambda_max / n,
          'ndebug': 0
        };
        (function() {
          var temp = runge_kutta.geodesic_simple(256, 1, x, v, opt);
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
        rf = final_x[1];
        delta_r = r0 - rf;
        m = 0.5; /* coordinates are such that mass=1/2 */
        g = ((m) * (Math.pow((r0), (-2.0)))); /* approximate as constant accel */
        delta_r_newtonian = ((0.5) * (g) * (Math.pow((lambda_max), (2.0))));
        rel_err = (delta_r - delta_r_newtonian) / delta_r_newtonian;
        expect_err = 2 * delta_r_newtonian / r0; /* expected rel error due to constant-accel approximation */
        if (verbosity >= 2) {
          print("delta_r=", delta_r, ", delta_r_newtonian=", delta_r_newtonian, " rel err=", rel_err, ", expected=", expect_err);
        }
        test.assert_rel_equal_eps(delta_r, delta_r_newtonian, expect_err);
        if (Math.abs(rel_err) > 2 * expect_err) {
          throw 'error in final r greater than expected';;
        }
        v_newtonian = -Math.sqrt(2 * m * (1 / rf - 1 / r0));
        vf = final_v[1];
        rel_err_v = (vf - v_newtonian) / v_newtonian;
        if (verbosity >= 2) {
          print("vf=", vf, ", vf_newtonian=", v_newtonian, " rel err=", rel_err_v, ", expected=", expect_err);
        }
        test.assert_rel_equal_eps(vf, v_newtonian, expect_err);
        if (Math.abs(rel_err_v) > 2 * expect_err) {
          throw 'error in final v greater than expected';;
        }
      };
      /*-------------------------------------------------------------------------------------------------- */
      test_runge_kutta.circular_orbit_period = function() {
        var r, v_phi, x, v, n, period, opt, err, final_x, final_v, final_a, final_lambda, info, eps;

        /*
        Period of a circular orbit, Schwarzschild coordinates.
        */
        r = 3.0;
        v_phi = 1 / Math.sqrt(2.0 * r * r * r); /* exact condition for circular orbit in Sch., if v_t=1. */
        x = [0.0, r, 1.0, 0.0, 0.0];
        v = [1.0, 0.0, 0.0, v_phi, 0.0];
        n = 100;
        period = 2.0 * Math.PI / v_phi;
        opt = {
          'lambda_max': period,
          'dlambda': period / n,
          'ndebug': 0,
          'norm_final': (false)
        };
        (function() {
          var temp = runge_kutta.geodesic_simple(256, 1, x, v, opt);
          err = temp[0];
          final_x = temp[1];
          final_v = temp[2];
          final_a = temp[3];
          final_lambda = temp[4];
          info = temp[5]
        })();
        if (verbosity >= 2) {
          print("final x=", io_util.vector_to_str_n_decimals(final_x, 16));
        }
        if (err & 1) {
          throw 'error: ' + info['message'];;
        }
        eps = ((10000.0) * (Math.pow((n), (-4.0))));
        test.assert_equal_eps(x[2], final_x[2], eps);
        test.assert_equal_eps(x[3], final_x[3], eps);
        test.assert_equal_eps(x[4], final_x[4], eps);
      };
      /*-------------------------------------------------------------------------------------------------- */
      test_runge_kutta.elliptical_orbit_period = function(r, a, direction, n, half_period) {
        var spacetime, chart, v_phi, x, circular_period, v, q, r_max, foo, period, triggers, lambda_max, ndebug, opt, err, final_x, final_v, final_a, final_lambda, info, eps, lamx, final_t, final_r, final_j;

        /*
        Period of an elliptical orbit, Schwarzschild coordinates.
        Start at perihelion, r. Make the initial velocity greater than the circular-orbit value by the factor a.
        Test against the Keplerian period. There is no point in testing with large n, because the errors
        become dominated by the Keplerian approximation.
        direction = angle about the x axis for the initial motion, defines plane of orbit
        Runge-Kutta with n steps.
        This is intended to be used with very large r, so that the Keplerian approximation is good.
        */
        spacetime = 256;
        chart = 1;
        v_phi = 1 / Math.sqrt(2.0 * r * r * r); /* exact condition for circular orbit in Sch., if v_t=1. */
        x = [0.0, r, 1.0, 0.0, 0.0];
        circular_period = 2.0 * Math.PI / v_phi;
        /* Increase velocity at perihelion to make orbit elliptical: */
        v_phi = v_phi * a;
        v = [1.0, 0.0, 0.0, v_phi * Math.cos(direction), v_phi * Math.sin(direction)];
        v = angular.make_tangent(x, v);
        v = vector.normalize(spacetime, chart, x, v);
        /* Compute newtonian r_max: */
        q = ((-1.0) + (((0.5) * (Math.pow((a), (2.0))))));
        r_max = ((((1.0) / (2.0))) * (Math.pow((q), (-1.0))) * (((-1.0) + (((-1.0) * (Math.pow((((1.0) + (((2.0) * (Math.pow((a), (2.0))) * (q))))), (((1.0) / (2.0))))))))) * (r));
        foo = ((1.0) + (((2.0) * (Math.pow((a), (2.0))) * (q))));
        period = ((0.35355339059327384) * (circular_period) * (Math.pow((((Math.pow((r), (-1.0))) * (((r) + (r_max))))), (1.5)))); /* Kepler's law of periods */
        triggers = [];
        lambda_max = period;
        if (half_period) {
          lambda_max = 0.5 * period * (1 + 2.0 / n); /* has to be longer than 0.5*period, or it doesn't test trigger */
          triggers = [
            [-1.0, 6, 0, 0.5]
          ]; /* trigger when rdot crosses zero from above */
        }
        /*-- */
        ndebug = 0;
        if (verbosity >= 3) {
          if (half_period) {
            print("testing half-period");
          }
          ndebug = n / 10;
        }
        opt = {
          'lambda_max': lambda_max,
          'dlambda': lambda_max / n,
          'ndebug': ndebug,
          'norm_final': (false),
          'triggers': triggers
        };
        (function() {
          var temp = runge_kutta.geodesic_simple(spacetime, chart, x, v, opt);
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
        if (verbosity >= 2) {
          print("final x=", io_util.vector_to_str_n_decimals(final_x, 16));
        }
        if (half_period) {
          eps = 10.0 / n; /* won't be as accurate as RK, because we extrapolate linearly to find where v crosses zero */
          lamx = -final_v[1] / final_a[1]; /* additional increment to lambda based on extrapolation to the trigger */
          final_lambda = final_lambda + lamx;
          final_t = final_x[0] + final_v[0] * lamx;
          final_r = final_x[1]; /* any correction would be second order */
          final_j = final_x[3] + final_v[3] * lamx;
          if (verbosity >= 3) {
            print("final lam=", final_lambda, ", t=", final_t, ", r=", final_r, ", j=", final_j);
          }
          test.assert_rel_equal_eps(final_r, r_max, eps);
          test.assert_rel_equal_eps(final_t, 0.5 * period, eps);
          test.assert_equal_eps(final_j, 0.0, eps);
        } else {
          eps = ((((10000.0) * (Math.pow((n), (-4.0))))) + (((100.0) * (Math.pow((r), (-1.0)))))); /* first term is for error in Keplerian period, second for Runge-Kutta */
          test.assert_equal_eps(x[2], final_x[2], eps);
          test.assert_equal_eps(x[3], final_x[3], eps);
          test.assert_equal_eps(x[4], final_x[4], eps);
        }
      };
      /*-------------------------------------------------------------------------------------------------- */
      /*verbosity=3 */
      test_runge_kutta.main = function() {
        var r, a, direction, n;

        test_runge_kutta.smoke_test();
        test_runge_kutta.simple_newtonian_free_fall();
        test_runge_kutta.circular_orbit_period();
        /*-- */
        r = 1.0e8;
        a = 1.1;
        direction = 0.0;
        n = 100;
        test_runge_kutta.elliptical_orbit_period(r, a, direction, n, (false)); /* test period */
        test_runge_kutta.elliptical_orbit_period(r, a, direction, n, (true)); /* test half-period */
        /*-- */
        test.done(verbosity, "runge_kutta");
      };
      test_runge_kutta.main();
