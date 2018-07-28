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
      /* ... exited due to a trigger */

      karl.load("schwarzschild");
      karl.load("kruskal");
      karl.load("keplerian");
      karl.load("angular");
      karl.load("transform");
      runge_kutta.trajectory_simple = function(spacetime, chart, pars, x0, v0, opt) {
        var x, v, lambda_max, dlambda, ndebug, debug_function, lambda0, norm_final, n_triggers, trigger_s, trigger_on, trigger_threshold, trigger_alpha, force_acts, force_function, user_function, force_chart, n, steps_between_debugging, debug_count, lam, ok, ndim, christoffel_function, use_c, ndim2, order, acc, y0, est, step, i, user_data, y, tr, tot_est;

        /*
        Calculate a trajectory using geodesic equation plus external force term, with 4th-order Runge-Kutta.
        spacetime = label for the spacetime we're doing (see spacetimes.h for labels)
        chart = label for the coordinate chart we're using
        x = starting point, given as an array of coordinates
        v = components of starting tangent vector
          The velocity need not be normalized, can be null or spacelike. If it is normalized, then the
          affine parameter represents proper time.
        The following options are in the hash opt[].
          lambda_max = maximum affine parameter, i.e., where to stop (but could stop earlier, e.g., if 
                         we hit a singularity)
          dlambda = step size
          ndebug = 0, or, if nonzero, determines how often to print debugging output; e.g., if ndebug=100
                     then we print debugging information at every 100th step
          debug_function = a function to be called when printing debugging output
          lambda0 = initial affine parameter, defaults to 0
          norm_final = adjust the final x and v to lie on and tangent to the unit sphere in i-j-k space;
                       default=(true)
          triggers = array of 4-element arrays, each describing a trigger (see below)
          force_acts = boolean, do we have an external force?
          force_function = function that calculates the proper acceleration vector d^2x/dlambda^2,
                                   given (lambda,x,v) as inputs; its output will be used and immediately discarded,
                                   so the function does not need to clone it before returning it
          force_chart = chart that the function wants for its inputs and outputs
          user_function = an optional user-defined function that is called on every iteration
        triggers
          These allow the integration to be halted when it appears that in the next iteration,
          a certain coordinate or velocity would cross a certain threshold.
          [0] = sense, +1 or -1 for triggering in a rising or falling direction
          [1] = index of coordinate (0-4) or velocity (5-9) on which to trigger
          [2] = threshold value
          [3] = fudge factor alpha, which should be less than 1 if something bad happens at threshold or if it would
                         be bad not to get the trigger
        returns
          [err,final_x,final_v,final_a,final_lambda,info]
        where
          err = 0 if normal, or bitwise or of codes such as 1, 2, defined in runge_kutta.h
          final_x,final_v,final_a,final_lambda = final values of position, velocity, acceleration, and affine param
          info = hash with misc. info, including error message and which trigger was triggered
        */
        x = (karl.clone_array1d(x0));
        v = (karl.clone_array1d(v0));
        /*-- process input options */
        (function() {
          var temp = runge_kutta.runge_kutta_get_options_helper(opt);
          lambda_max = temp[0];
          dlambda = temp[1];
          ndebug = temp[2];
          debug_function = temp[3];
          lambda0 = temp[4];
          norm_final = temp[5];
          n_triggers = temp[6];
          trigger_s = temp[7];
          trigger_on = temp[8];
          trigger_threshold = temp[9];
          trigger_alpha = temp[10];
          force_acts = temp[11];
          force_function = temp[12];
          user_function = temp[13];
          force_chart = temp[14]
        })();
        /*-- initial setup */
        (function() {
          var temp = runge_kutta.runge_kutta_init_helper(lambda_max, lambda0, dlambda, ndebug, spacetime, chart, pars);
          n = temp[0];
          steps_between_debugging = temp[1];
          debug_count = temp[2];
          lam = temp[3];
          ok = temp[4];
          ndim = temp[5];
          christoffel_function = temp[6]
        })();
        if (!ok) {
          return [1, x, v, 0.0, runge_kutta.mess(["unrecognized spacetime or chart: ", spacetime, " ", chart])];
        }
        use_c = runge_kutta.c_available(spacetime, chart);
        ndim2 = ndim * 2; /* Reduce 2nd-order ODE to ndim2 coupled 1st-order ODEs. */
        if (((x).length) != ndim || ((v).length) != ndim) {
          return [1, x, v, 0.0, runge_kutta.mess(["x or v has wrong length"])];
        }
        order = 4; /* 4th order Runge-Kutta */
        acc = karl.array1d((ndim));
        y0 = karl.array1d((ndim2));
        if (use_c) {
          /* use faster C implementation: */
        } else {
          runge_kutta.find_acc = function(acc, y, dlambda) {
            runge_kutta.apply_christoffel(christoffel_function, y, acc, dlambda, ndim);
          }
        }
        runge_kutta.next_est = function(est, acc, step, y, dlambda) {
          runge_kutta.find_acc(acc, y, dlambda);
          if (force_acts) {
            runge_kutta.handle_force(acc, lam, x, v, force_function, force_chart, ndim, spacetime, chart, pars, dlambda);
          }
          for (var i = 0; i < ndim; i++) {
            est[step][i] = y[ndim + i] * dlambda;
            est[step][ndim + i] = acc[i];
          }
        }
        user_data = null;
        for (var iter = 0; iter < n; iter++) {
          dlambda = (lambda_max - lam) / (n - iter); /* small readjustment so we land on the right final lambda */
          est = karl.array2d(ndim2, order);
          /*         =k in the notation of most authors */
          /*         Four estimates of the changes in the independent variables for 4th-order Runge-Kutta. */
          debug_count = runge_kutta.debug_helper(debug_count, ndebug, steps_between_debugging, iter, lam, dlambda, x, v, debug_function, spacetime | chart);
          for (var i = 0; i < ndim; i++) {
            y0[i] = x[i];
          }
          for (var i = 0; i < ndim; i++) {
            y0[i + ndim] = v[i];
          }
          y0 = (karl.clone_array1d(y0));
          /* ...Disentangle it from x and v so that in python, changing x or v can't change it. This is actually */
          /*    not necessary, because x and v don't change until the next iteration, when y0 is built again, */
          /*    but I find it too hard to reason about the code without this. */
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
            runge_kutta.next_est(est, acc, step, y, dlambda);
          }
          if (n_triggers > 0) {
            tr = runge_kutta.trigger_helper(x, v, acc, dlambda, n_triggers, trigger_s, trigger_on, trigger_threshold, trigger_alpha, ndim);
            if (tr != (-1)) {
              return runge_kutta.runge_kutta_final_helper(debug_count, ndebug, steps_between_debugging, iter, lam, dlambda, x, v, acc, norm_final, debug_function, spacetime | chart, pars, user_data, 3, tr);
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
          /*-- Call user function: */
          if (!((user_function) == null)) {
            user_data = user_function(lam, x, v, spacetime, chart, pars);
          }
        }
        return runge_kutta.runge_kutta_final_helper(debug_count, ndebug, steps_between_debugging, n, lam, dlambda, x, v, acc, norm_final, debug_function, spacetime | chart, pars, user_data, 0, -1);
      };
      runge_kutta.c_available = function(spacetime, chart) {
        return false;
      };
      runge_kutta.unpack_pars = function(pars_array, spacetime, chart) {

        /*
        I'm trying to keep the interface to the C code simple, so for its use, the parameters in the pars
        hash are unpacked into an array of floats. The array of floats is allocated by the calling code
        as a C array with a fixed size.
        */
        if (spacetime == 256) {
          return; /* Schwarzschild spacetime has no adjustable parameters. */
        }
        throw "unimplemented spacetime: " + str(spacetime);;
      };
      runge_kutta.r_stuff = function(spacetime, chart, pars, x, v, acc, pt, acc_p, pt_p) {
        var ndim, ndim2, x2, r, v2, pt, i, ok, christoffel_function, name, rdot, rddot, p, lam_left;

        /*
        Returns [err,r,r',r'',p,lam_left], where the primes represent derivatives with respect to affine parameter,
        and for small r,
        p is an estimate of the exponent in r ~ lambda^p, and lam_left is an estimate of the distance
        left before the singularity in terms of the affine parameter. For values of r that are not small,
        the p and lam_left outputs would be meaningless, and are returned as NaN.
        The arrays acc_p and pt_p are pointers to arrays that have already been allocated,
        or None if this is the js implementation. If r<0, err=1.
        If r==1, then we return [0,r,(NaN),(NaN),(NaN),(NaN)].
        */
        ndim = 5;
        ndim2 = 10;
        x2 = transform.transform_point(x, spacetime, chart, pars, 1);
        r = x2[1];
        if ((isNaN(r)) || r < 0.0) {
          return [1, r, 0.0, 0.0, 0.0, 0.0];
        }
        if (r == 1.0) {
          /* We'd get an error below in transforming the vector to Schwarzschild. */
          return [0, r, (NaN), (NaN), (NaN), (NaN)];
        }
        v2 = transform.transform_vector(v, x, spacetime, chart, pars, 1);
        for (var i = 0; i < ndim; i++) {
          pt[i] = x2[i];
          pt[i + 5] = v2[i];
        }
        if (runge_kutta.c_available(256, 1)) {
          /* use faster C implementation: */
        } else {
          (function() {
            var temp = transform.chart_info(256 | 1);
            ok = temp[0];
            ndim = temp[1];
            christoffel_function = temp[2];
            name = temp[3]
          })();
          runge_kutta.apply_christoffel(christoffel_function, pt, acc, 1.0, ndim);
        }
        rdot = v2[1];
        rddot = acc[1];
        p = (NaN);
        lam_left = (NaN);
        if (r < 1.0) {
          p = 1.0 / (1 - r * rddot / (rdot * rdot));
          if (p < 0.4) {
            p = 0.4;
          }
          if (p > 1.0) {
            p = 1.0;
          }
          lam_left = -p * r / rdot; /* estimate of when we'd hit the singularity */
        }
        return [0, r, rdot, rddot, p, lam_left];
      };
      runge_kutta.handle_force = function(a, lam, x, v, force_function, force_chart, ndim, spacetime, chart, pars, dlambda) {
        var x2, v2, proper_accel2, proper_accel, a, i;

        /* The API says that force_function does not need to clone its output vector before returning it, so */
        /* we need to make sure to discard it here and never do anything with it later. */
        /* The vector a that we're modifying already has a factor of dlambda in it, so we multiply by dlambda here */
        /* as well. */
        x2 = transform.transform_point(x, spacetime, chart, pars, force_chart);
        v2 = transform.transform_vector(v, x, spacetime, chart, pars, force_chart);
        proper_accel2 = force_function(lam, x2, v2);
        proper_accel = transform.transform_vector(proper_accel2, x2, spacetime, force_chart, pars, chart);
        for (var i = 0; i < ndim; i++) {
          a[i] = a[i] + proper_accel[i] * dlambda;
        }
      };
      runge_kutta.runge_kutta_get_options_helper = function(opt) {
        var lambda_max, dlambda, ndebug, debug_function, lambda0, norm_final, n_triggers, trigger_s, trigger_on, trigger_threshold, trigger_alpha, user_function, force_acts, force_function, force_chart;

        lambda_max = runge_kutta.runge_kutta_get_par_helper(opt, "lambda_max", null);
        dlambda = runge_kutta.runge_kutta_get_par_helper(opt, "dlambda", null);
        ndebug = runge_kutta.runge_kutta_get_par_helper(opt, "ndebug", 0);
        debug_function = runge_kutta.runge_kutta_get_par_helper(opt, "debug_function", runge_kutta.default_debug_function);
        lambda0 = runge_kutta.runge_kutta_get_par_helper(opt, "lambda0", 0.0);
        norm_final = runge_kutta.runge_kutta_get_par_helper(opt, "norm_final", (true));
        n_triggers = 0;
        (function() {
          var temp = [
            [],
            [],
            [],
            []
          ];
          trigger_s = temp[0];
          trigger_on = temp[1];
          trigger_threshold = temp[2];
          trigger_alpha = temp[3]
        })();
        if ((("triggers") in (opt))) {
          n_triggers = runge_kutta.runge_kutta_get_trigger_options_helper(opt, trigger_s, trigger_on, trigger_threshold, trigger_alpha);
        }
        user_function = null;
        if ((("user_function") in (opt))) {
          user_function = opt['user_function'];
        }
        force_acts = runge_kutta.runge_kutta_get_par_helper(opt, "force_acts", (false));
        force_function = runge_kutta.runge_kutta_get_par_helper(opt, "force_function", 0);
        force_chart = runge_kutta.runge_kutta_get_par_helper(opt, "force_chart", 0);
        return [lambda_max, dlambda, ndebug, debug_function, lambda0, norm_final, n_triggers, trigger_s, trigger_on, trigger_threshold, trigger_alpha, force_acts, force_function, user_function, force_chart];
      };
      runge_kutta.runge_kutta_init_helper = function(lambda_max, lambda0, dlambda, ndebug, spacetime, chart, pars) {
        var n, steps_between_debugging, debug_count, lam, ok, ndim, christoffel_function, name;

        n = Math.ceil((lambda_max - lambda0) / dlambda); /* dlambda will be adjusted slightly in order to deal with the rounding */
        if (ndebug == 0) {
          steps_between_debugging = n * 2; /* debugging will never happen */
        } else {
          steps_between_debugging = ndebug;
        }
        debug_count = steps_between_debugging + 1; /* trigger it on the first iteration */
        lam = lambda0;
        (function() {
          var temp = transform.chart_info(spacetime | chart);
          ok = temp[0];
          ndim = temp[1];
          christoffel_function = temp[2];
          name = temp[3]
        })();
        return [n, steps_between_debugging, debug_count, lam, ok, ndim, christoffel_function];
      };
      runge_kutta.trigger_helper = function(x, v, acc, dlambda, n_triggers, trigger_s, trigger_on, trigger_threshold, trigger_alpha, ndim) {
        var s, m, thr, alpha, dx, x_dot;

        /*
        Check triggers:
        We can trigger in the rising (s=+1) or falling (s=-1) direction. The coordinate or velocity
        we're triggering on differs from the trigger value by dx, and it's currently changing at
        a rate x_dot. Depending on the signs of s, dx, and x_dot, we have 8 cases. The logic below
        handles all the cases properly. The input variable acc is actually the acceleration times dlambda.
        Returns -1 normally, or trigger number if a trigger tripped off.
        */
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
            x_dot = acc[m - ndim] / dlambda; /* left over from step==3, good enough for an estimate */
          }
          /*print("s=",s,", dx=",dx,", x_dot=",x_dot,", dlambda=",dlambda,"lhs=",x_dot*dlambda*s,", rhs=",alpha*dx*s) */
          if (s * dx > 0 && x_dot * dlambda * s > alpha * dx * s) { /* Note that we can't cancel the s, which may be negative. */
            /* We extrapolate that if we were to complete this iteration, we would cross the threshold. */
            return i;
          }
        }
        return -1;
      };
      runge_kutta.runge_kutta_final_helper = function(debug_count, ndebug, steps_between_debugging, n, lam, dlambda, x, v, acc, norm_final, debug_function, spacetime_and_chart, pars, user_data, err, which_trigger) {
        var x, v, info;

        runge_kutta.debug_helper(debug_count, ndebug, steps_between_debugging, n, lam, dlambda, x, v, debug_function, spacetime_and_chart);
        /* ... always do a printout for the final iteratation */
        if (norm_final) {
          x = angular.renormalize(x);
          v = angular.make_tangent(x, v);
        }
        info = {};
        if (!((user_data) == null)) {
          info['user_data'] = user_data;
        }
        if (which_trigger != -1) {
          info['trigger'] = which_trigger;
        }
        return [err, x, v, acc, lam, info];
      };
      runge_kutta.runge_kutta_get_par_helper = function(opt, name, default_val) {
        var val;

        val = default_val;
        if (((name) in (opt))) {
          val = opt[name];
        }
        if (((val) == null)) {
          throw 'required option ' + name + '  !  supplied';;
        }
        return val;
      };
      runge_kutta.runge_kutta_get_trigger_options_helper = function(opt, trigger_s, trigger_on, trigger_threshold, trigger_alpha) {
        var n_triggers, trigger;

        n_triggers = ((opt["triggers"]).length);
        for (var i = 0; i < n_triggers; i++) {
          trigger = opt["triggers"][i];
          ((trigger_s).push(trigger[0]));
          ((trigger_on).push(trigger[1]));
          ((trigger_threshold).push(trigger[2]));
          ((trigger_alpha).push(trigger[3]));
        }
        return n_triggers;
      };
      runge_kutta.apply_christoffel = function(christoffel_function, y, acc, dlambda, ndim) {
        var ch, a, acc, i;

        /* A similar routine, written in C for speed, is in apply_christoffel.c. This version exists */
        /* only so it can be translated into javascript. */
        ch = christoffel_function(y);
        for (var i = 0; i < ndim; i++) {
          a = 0.0; /* is essentially the acceleration */
          for (var j = 0; j < ndim; j++) {
            for (var k = 0; k < ndim; k++) {
              a -= ch[j][k][i] * y[ndim + j] * y[ndim + k];
            }
          }
          acc[i] = a * dlambda;
        }
      };
      runge_kutta.mess = function(stuff) {

        return {
          'message': io_util.strcat(stuff)
        };
      };
      runge_kutta.debug_helper = function(debug_count, ndebug, steps_between_debugging, iter, lam, dlambda, x, v, debug_function, spacetime_and_chart) {
        var do_debug, debug_count, ok, ndim, christoffel_function, name;

        /*
        Prints debugging info and returns the updated debug_count.
        */
        do_debug = (false);
        if (ndebug != 0 && (debug_count >= steps_between_debugging)) {
          debug_count = 0;
          do_debug = (true);
        }
        if (do_debug) {
          (function() {
            var temp = transform.chart_info(spacetime_and_chart) /* just to get name of chart */ ;
            ok = temp[0];
            ndim = temp[1];
            christoffel_function = temp[2];
            name = temp[3]
          })();
          debug_function(iter, lam, dlambda, x, v, name);
        }
        return debug_count + 1;
      };
      runge_kutta.default_debug_function = function(iter, lam, dlambda, x, v, name) {

        print("i=", iter, " lam=", io_util.fl(lam), " x=", io_util.vector_to_str_n_decimals(x, 1), " v=", io_util.vector_to_str_n_decimals(v, 1));;
      };
