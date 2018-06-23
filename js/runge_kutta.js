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
      /* ... sch5 coordinates */
      /* ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...),  */
      /* return codes for Runge-Kutta, designed to be bitwise or-able. */
      /* ... something went really wrong, output is garbage */
      /* ... the geodesic was incomplete */

      karl.load("schwarzschild");
      karl.load("kruskal");
      karl.load("angular");
      runge_kutta.geodesic_simple = function(spacetime, chart, x0, v0, opt) {
          var x, v, lambda_max, dlambda, ndebug, lambda0, norm_final, n, do_limit_change, limit_change, ok, steps_between_debugging, debug_count, lam, ndim, christoffel_function, ndim2, order, y0, i, y, ch, a, est, step, tot_est;

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
            do_limit_change = boolean, should we do a sanity check by limiting changes in coordinates per step,
                    die if limit is violated?; default=(false)
            limit_change = this is approximately the maximum fractional change in r or 1/10 the maximum change in i,j,k,
                    expressed in units of 1/n; default: 1
          returns
            [err,final_x,final_v,final_lambda,info]
          where
            err = 0 if normal, or bitwise or of codes such as 1, 2, defined in runge_kutta.h
            final_x,final_v,final_lambda = final values from integration
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
          n = Math.ceil(lambda_max / dlambda);
          do_limit_change = (false);
          if ((("do_limit_change") in (opt))) {
              do_limit_change = opt["do_limit_change"];
          }
          if (do_limit_change) {
              limit_change = 1.0 / n;
              if ((("limit_change") in (opt))) {
                  do_limit_change = float(opt["do_limit_change"]) / n;
              }
          }
          ok = (false);
          if (ndebug == 0) {
              steps_between_debugging = n * 2; /* debugging will never happen */
          } else {
              steps_between_debugging = ndebug;
          }
          debug_count = steps_between_debugging + 1; /* trigger it on the first iteration */
          lam = 0.0;
          (function() {
              var temp = runge_kutta.chart_info(spacetime, chart);
              ok = temp[0];
              ndim = temp[1];
              christoffel_function = temp[2]
          })();
          if (!ok) {
              return [1, x, v, 0.0, runge_kutta.mess(["unrecognized spacetime or chart: ", spacetime, " ", chart])];
          }
          ndim2 = ndim * 2; /* Reduce 2nd-order ODE to ndim2 coupled 1st-order ODEs. */
          if (((x).length) != ndim || ((v).length) != ndim) {
              return [1, x, v, 0.0, runge_kutta.mess(["x or v has wrong length"])];
          }
          order = 4; /* 4th order Runge-Kutta */
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
                  ch = christoffel_function(y);
                  for (var i = 0; i < ndim2; i++);
                  for (var i = 0; i < ndim; i++) {
                      a = 0.0; /* is essentially the acceleration */
                      for (var j = 0; j < ndim; j++) {
                          for (var k = 0; k < ndim; k++) {
                              a -= ch[j][k][i] * y[ndim + j] * y[ndim + k];
                          }
                      }
                      est[step][ndim + i] = a * dlambda;
                  }
                  for (var i = 0; i < ndim; i++) {
                      est[step][i] = y[ndim + i] * dlambda;
                  }
              }
              lam = lam + dlambda;
              tot_est = karl.array1d((ndim2));
              for (var i = 0; i < ndim2; i++) {
                  tot_est[i] = (est[0][i] + 2.0 * est[1][i] + 2.0 * est[2][i] + est[3][i]) / 6.0;
              }
              for (var i = 0; i < ndim; i++) {
                  v[i] += tot_est[ndim + i];
              }
              for (var i = 0; i < ndim; i++) {
                  if (do_limit_change) {
                      runge_kutta.check_limit_change(spacetime, chart, x, tot_est, limit_change);
                  }
                  x[i] += tot_est[i];
              }
          }
          runge_kutta.debug_helper(debug_count, ndebug, steps_between_debugging, n, lam, x, v);
          if (norm_final) {
              x = angular.renormalize(x);
              v = angular.make_tangent(x, v);
          }
          return [0, x, v, lam, {}];
      };
      runge_kutta.check_limit_change = function(spacetime, chart, x, dx, limit_change) {
          var ok, rel_dr;

          /*
          Sanity check to flag sudden large changes in coordinates.
          */
          ok = (true);
          if ((spacetime | chart) == (256 | 1)) {
              rel_dr = Math.abs(dx[1]) / x[1];
          }
          if ((spacetime | chart) == (256 | 2)) {
              rel_dr = (Math.abs(dx[0]) + Math.abs(dx[1])) / (1 + Math.abs(x[0] - x[1]));
              /* ... quick and dirty estimate using r=a-b+1, not really appropriate for small distances */
          }
          if (rel_dr > limit_change) {
              throw io_util.strcat(['r changed by too much , rel_dr=', rel_dr, ', x=', io_util.vector_to_str(x), ', dx=', io_util.vector_to_str(dx)]);;
          }
          for (var i = 2; i < 5; i++) {
              if (Math.abs(dx[i]) > 10.0 * limit_change) {
                  throw 'angular coord. changed by too much';;
              }
          }
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
              print("i=", iter, " lam=", ("%4.2e" % lam), " x=", io_util.vector_to_str_n_decimals(x, 1), " v=", io_util.vector_to_str_n_decimals(v, 1));
          }
          return debug_count + 1;;
      };
