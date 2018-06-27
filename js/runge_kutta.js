      /*
               --- module runge_kutta ---
               This was translated from python. Do not edit directly.
            */
      if (typeof runge_kutta === 'undefined') {
        var runge_kutta = {};
      }

      ;
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

      karl.load("schwarzschild");
      karl.load("kruskal");
      karl.load("angular");
      runge_kutta.geodesic_simple = function(spacetime, chart, x0, v0, opt) {
        var x, v, lambda_max, dlambda, ndebug, lambda0, norm_final, n, ok, steps_between_debugging, n_triggers, trigger_s, trigger_on, trigger_threshold, trigger_alpha, trigger, i, debug_count, lam, ndim, christoffel_function, ndim2, order, acc, y0, y, ch, a, est, step, s, m, thr, alpha, dx, x_dot, tot_est;

        /*
        Calculate a geodesic using geodesic equation and 4th-order Runge-Kutta.
        spacetime = label for the spacetime we're doing (see spacetimes.h for labels)
        chart = label for the coordinate chart we're using
        x = starting point, given as an array of coordinates
        v = components of starting tangent vector (need not be normalized, can be null or spacelike)
        The following options are in the hash opt[].
          lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
                         we hit a singularity)
          dlambda = step size
          ndebug = 0, or, if nonzero, determines how often to print debugging output; e.g., if ndebug=100
                     then we print debugging information at every 100th step
          lambda0 = initial affine parameter, defaults to 0
          norm_final = adjust the final x and v to lie on and tangent to the unit sphere in i-j-k space;
                       default=(true)
          triggers = array of 4-element arrays, each describing a trigger (see below)
        triggers
          These allow the integration to be halted when it appears that in the next iteration,
          a certain coordinate or velocity would cross a certain threshold.
          [0] = sense, +1 or -1 for triggering in a rising or falling direction
          [1] = index of coordinate (0-4) or velocity (5-9) on which to trigger
          [2] = threshold value
          [3] = fudge factor alpha, which can should typically be less than 1; see docs
        returns
          [err,final_x,final_v,final_a,final_lambda,info]
        where
          err = 0 if normal, or bitwise or of codes such as 1, 2, defined in runge_kutta.h
          final_x,final_v,final_a,final_lambda = final values of position, velocity, acceleration, and affine param
          info = hash with keys below
            message = error message
        */
        x = (karl.clone_array1d(x0));
        v = (karl.clone_array1d(v0));
        lambda_max = opt["lambda_max"];
        dlambda = opt["dlambda"];
        ndebug = opt["ndebug"];
        lambda0 = 0.0;
        if ((("lambda0") in (opt))) {
          lambda0 = opt["lambda0"];
        }
        norm_final = (true);
        if ((("norm_final") in (opt))) {
          norm_final = opt["norm_final"];
        }
        n = Math.ceil((lambda_max - lambda0) / dlambda);
        ok = (false);
        if (ndebug == 0) {
          steps_between_debugging = n * 2; /* debugging will never happen */
        } else {
          steps_between_debugging = ndebug;
        }
        n_triggers = 0;
        if ((("triggers") in (opt))) {
          n_triggers = ((opt["triggers"]).length);
          trigger_s = [];
          trigger_on = [];
          trigger_threshold = [];
          trigger_alpha = [];
          for (var i = 0; i < n_triggers; i++) {
            trigger = opt["triggers"][i];
            trigger_s[i] = [0];
            trigger_on[i] = [1];
            trigger_threshold[i] = [2];
            trigger_alpha[i] = [3];
          }
        }
        debug_count = steps_between_debugging + 1; /* trigger it on the first iteration */
        lam = lambda0;
        (function() {
          var temp = runge_kutta.chart_info(spacetime, chart);
          ok = temp[0];
          ndim = temp[1];
          christoffel_function = temp[2]
        })();
        use_c = false;
        if (!ok) {
          return [1, x, v, 0.0, runge_kutta.mess(["unrecognized spacetime or chart: ", spacetime, " ", chart])];
        }
        ndim2 = ndim * 2; /* Reduce 2nd-order ODE to ndim2 coupled 1st-order ODEs. */
        if (((x).length) != ndim || ((v).length) != ndim) {
          return [1, x, v, 0.0, runge_kutta.mess(["x or v has wrong length"])];
        }
        order = 4; /* 4th order Runge-Kutta */
        acc = karl.array1d((ndim));
        for (var iter = 0; iter < n; iter++) {
          est = karl.array2d(ndim2, order);; /*         =k in the notation of most authors */
          /*         Four estimates of the changes in the independent variables for 4th-order Runge-Kutta. */
          debug_count = runge_kutta.debug_helper(debug_count, ndebug, steps_between_debugging, iter, lam, x, v);
          y0 = karl.array1d((ndim2));
          for (var i = 0; i < ndim; i++) {
            y0[i] = (x[i]);
          }
          for (var i = 0; i < ndim; i++) {
            y0[i + ndim] = (v[i]);
          }
          for (var step = 0; step < order; step++) {
            if (step == 0) {
              y = (karl.clone_array1d(y0));
            }
            if (step == 1) {
              for (var i = 0; i < ndim2; i++) {
                y[i] = y0[i] + 0.5 * est[0][i];
              }
            }
            if (step == 2) {
              for (var i = 0; i < ndim2; i++) {
                y[i] = y0[i] + 0.5 * est[1][i];
              }
            }
            if (step == 3) {
              for (var i = 0; i < ndim2; i++) {
                y[i] = y0[i] + est[2][i];
              }
            }
            for (var i = 0; i < ndim2; i++);
            if (use_c) {
              /* use faster C implementation: */
            } else {
              ch = christoffel_function(y);
              for (var i = 0; i < ndim; i++) {
                a = 0.0; /* is essentially the acceleration */
                for (var j = 0; j < ndim; j++) {
                  for (var k = 0; k < ndim; k++) {
                    a -= ch[j][k][i] * y[ndim + j] * y[ndim + k];
                  }
                }
                est[step][ndim + i] = a * dlambda;
                acc[i] = a; /* may be needed for trigger detection */
              }
            }
            for (var i = 0; i < ndim; i++) {
              est[step][i] = y[ndim + i] * dlambda;
            }
          }
          /*-- Check triggers: */
          /* We can trigger in the rising (s=+1) or falling (s=-1) direction. The coordinate or velocity */
          /* we're triggering differs from the trigger value by dx, and it's currently changing at */
          /* a rate x_dot. Depending on the signs of s, dx, and x_dot, we have 8 cases. The logic below */
          /* handles all the cases properly. */
          for (var i = 0; i < n_triggers; i++) {
            s = trigger_s[i]; /* sense of the trigger (see above) */
            m = trigger_on[i]; /* index of coordinate or velocity on which to trigger */
            thr = trigger_threshold[i]; /* threshold value */
            alpha = trigger_alpha[i]; /* fudge factor, can typically be 1; see docs */
            if (m < ndim) {
              /* triggering on a coordinate */
              dx = thr - x[m];
              x_dot = v[m];
            } else {
              /* triggering on a velocity */
              dx = thr - v[m - ndim];
              x_dot = acc[m - ndim]; /* left over from step==3, good enough for an estimate */
            }
            if (s * dx > 0 && x_dot * dlambda * s > alpha * dx * s) { /* Note that we can't cancel the s, which may be negative. */
              /* We extrapolate that if we were to complete this iteration, we would cross the threshold. */
              return runge_kutta.runge_kutta_final_helper(debug_count, ndebug, steps_between_debugging, n, lam, x, v, acc, norm_final);
            }
          }
          /*-- Update everything: */
          lam = lam + dlambda;
          tot_est = karl.array1d((ndim2));
          for (var i = 0; i < ndim2; i++) {
            tot_est[i] = (est[0][i] + 2.0 * est[1][i] + 2.0 * est[2][i] + est[3][i]) / 6.0;
          }
          for (var i = 0; i < ndim; i++) {
            v[i] += tot_est[ndim + i];
          }
          for (var i = 0; i < ndim; i++) {
            x[i] += tot_est[i];
          }
        }
        return runge_kutta.runge_kutta_final_helper(debug_count, ndebug, steps_between_debugging, n, lam, x, v, acc, norm_final);
      };
      runge_kutta.runge_kutta_final_helper = function(debug_count, ndebug, steps_between_debugging, n, lam, x, v, acc, norm_final) {
        var x, v;

        runge_kutta.debug_helper(debug_count, ndebug, steps_between_debugging, n, lam, x, v);
        /* ... always do a printout for the final iteratation */
        if (norm_final) {
          x = angular.renormalize(x);
          v = angular.make_tangent(x, v);
        }
        return [0, x, v, acc, lam, {}];
      };
      runge_kutta.mess = function(stuff) {

        return {
          'runge_kutta.message': io_util.strcat(stuff)
        };
      };
      runge_kutta.chart_info = function(spacetime, chart) {
        var recognized;

        recognized = (false);
        if ((spacetime | chart) == (256 | 1)) {
          return [(true), 5, schwarzschild.christoffel];
        }
        if ((spacetime | chart) == (256 | 2)) {
          return [(true), 5, kruskal.christoffel];
        }
        return [(false), None, None];
      };
      runge_kutta.debug_helper = function(debug_count, ndebug, steps_between_debugging, iter, lam, x, v) {
        var do_debug, debug_count;

        /*
        Prints debugging info and returns the updated debug_count.
        */
        do_debug = (false);
        if (ndebug != 0 && (debug_count >= steps_between_debugging)) {
          debug_count = 0;
          do_debug = (true);
        }
        if (do_debug) {
          print("i=", iter, " lam=", io_util.fl(lam), " x=", io_util.vector_to_str_n_decimals(x, 1), " v=", io_util.vector_to_str_n_decimals(v, 1));
        }
        return debug_count + 1;;
      };
