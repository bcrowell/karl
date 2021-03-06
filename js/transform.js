      /*
               --- module transform ---
               This was translated from python. Do not edit directly.
            */
      if (typeof transform === 'undefined') {
        var transform = {};
      }

      /*
      Routines to compute transformations of points among Schwarzschild,
      arcsinh-Kruskal, and Keplerian coordinates, the jacobians of those transformations,
      and transformations of vectors.
      Documentation for the math is in the file doc.tex, which can be
      compiled to pdf format by doing a "make doc." (Comments in the code do
      not document the math or the definitions of the variables.) 
      */
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
      karl.load("math_util");
      karl.load("kruskal");
      karl.load("keplerian");
      karl.load("schwarzschild");
      karl.load("io_util");
      transform.chart_info = function(spacetime_and_chart) {
        var recognized;

        /*
        Input is spacetime|chart. Returns [ok,ndim,christoffel_function,name].
        See README for a list of things to do when adding a new chart.
        */
        recognized = (false);
        if ((spacetime_and_chart) == (256 | 1)) {
          return [(true), 5, schwarzschild.christoffel, "SCH"];
        }
        if ((spacetime_and_chart) == (256 | 2)) {
          return [(true), 5, kruskal.christoffel, "AKS"];
        }
        if ((spacetime_and_chart) == (256 | 3)) {
          return [(true), 5, keplerian.christoffel, "KEP"];
        }
        return [(false), None, None, ''];
      };
      transform.transform_point = function(x, spacetime, chart, pars, chart2) {
        var ok, ndim, christoffel_function, name, ndim2, x2, i, done, a, b, t, r;

        /*
        Transforms a point x from chart to chart2. Return value is an array, which is automatically cloned.
        It's legal to have chart==chart2. If the input point is not in the domain of the relevant functions,
        returns an array whose elements are all NaN.
        */
        if (chart == chart2) {
          return (karl.clone_array1d(x));
        }
        if (spacetime != 256) {
          throw io_util.strcat((["unrecognized spacetime=", spacetime]));;
        }
        (function() {
          var temp = transform.chart_info(spacetime | chart);
          ok = temp[0];
          ndim = temp[1];
          christoffel_function = temp[2];
          name = temp[3]
        })();
        if (!ok) {
          throw io_util.strcat((["unrecognized chart, spacetime=", spacetime, "chart=", chart]));;
        }
        (function() {
          var temp = transform.chart_info(spacetime | chart2);
          ok = temp[0];
          ndim2 = temp[1];
          christoffel_function = temp[2];
          name = temp[3]
        })();
        if (!ok) {
          throw io_util.strcat((["unrecognized chart, spacetime=", spacetime, "chart=", chart2]));;
        }
        x2 = karl.array1d((ndim2));
        if ((chart == 1 || chart == 3) && x[1] < 0.0) {
          for (var i = 0; i < ndim2; i++) {
            x2[i] = (NaN);
          }
          return (karl.clone_array1d(x2)); /* cloning not really needed in this case, but be consistent */
        }
        done = (false);
        /* The following assumes, as is true for 1, 2, and 3, that coords 2,3,4 are (i,j,k). */
        x2[2] = x[2];
        x2[3] = x[3];
        x2[4] = x[4];
        /* Transform the first two coordinates: */
        if (chart == 1 && chart2 == 2) {
          (function() {
            var temp = transform.schwarzschild_to_kruskal(x[0], x[1]);
            a = temp[0];
            b = temp[1]
          })();
          x2[0] = a;
          x2[1] = b;
          done = (true);
        }
        if (!done && chart == 2 && chart2 == 1) {
          (function() {
            var temp = transform.kruskal_to_schwarzschild(x[0], x[1]);
            t = temp[0];
            r = temp[1]
          })();
          x2[0] = t;
          x2[1] = r;
          done = (true);
        }
        if (!done && chart == 1 && chart2 == 3) {
          x2[0] = x[0];
          x2[1] = Math.pow((x[1]), (1.5));
          done = (true);
        }
        if (!done && chart == 3 && chart2 == 1) {
          x2[0] = x[0];
          x2[1] = Math.pow((x[1]), (0.6666666666666666));
          done = (true);
        }
        if (!done && chart == 2 && chart2 == 3) {
          (function() {
            var temp = transform.kruskal_to_schwarzschild(x[0], x[1]);
            t = temp[0];
            r = temp[1]
          })();
          x2[0] = t;
          x2[1] = Math.pow((r), (1.5));
          done = (true);
        }
        if (!done && chart == 3 && chart2 == 2) {
          t = x[0];
          r = Math.pow((x[1]), (0.6666666666666666));
          (function() {
            var temp = transform.schwarzschild_to_kruskal(t, r);
            a = temp[0];
            b = temp[1]
          })();
          x2[0] = a;
          x2[1] = b;
          done = (true);
        }
        if (done) {
          return (karl.clone_array1d(x2)); /* for python, divorce the new vector from entanglement with components of old */
        } else {
          throw io_util.strcat((["don't know how to transform, spacetime=", spacetime, ", chart=", chart, ", chart2=", chart2]));;
        }
      };
      transform.transform_vector = function(v, x, spacetime, chart, pars, chart2) {
        var ok, ndim, christoffel_function, name, ndim2, v2, found_jac, t, r, jac, u;

        /*
        Transforms a vector v from chart to chart2 in the tangent space at x (x being described in chart, not chart2).
        Return value is an array, which is
        automatically cloned. It's legal to have chart==chart2.
        */
        if (chart == chart2) {
          return (karl.clone_array1d(v));
        }
        if (spacetime != 256) {
          throw io_util.strcat((["unrecognized spacetime=", spacetime]));;
        }
        (function() {
          var temp = transform.chart_info(spacetime | chart);
          ok = temp[0];
          ndim = temp[1];
          christoffel_function = temp[2];
          name = temp[3]
        })();
        if (!ok) {
          throw io_util.strcat((["unrecognized chart, spacetime=", spacetime, "chart=", chart]));;
        }
        (function() {
          var temp = transform.chart_info(spacetime | chart2);
          ok = temp[0];
          ndim2 = temp[1];
          christoffel_function = temp[2];
          name = temp[3]
        })();
        if (!ok) {
          throw io_util.strcat((["unrecognized chart, spacetime=", spacetime, "chart=", chart2]));;
        }
        v2 = karl.array1d((ndim2));
        /* The following assumes, as is true for 1, 2, and 3, that coords 2,3,4 are (i,j,k). */
        v2[2] = v[2];
        v2[3] = v[3];
        v2[4] = v[4];
        found_jac = (false);
        /* Find the Jacobian. */
        if (chart == 1 && chart2 == 2) {
          (function() {
            var temp = [x[0], x[1]];
            t = temp[0];
            r = temp[1]
          })();
          jac = transform.jacobian_schwarzschild_to_kruskal(t, r);
          found_jac = (true);
        }
        if (chart == 2 && chart2 == 1) {
          (function() {
            var temp = transform.kruskal_to_schwarzschild(x[0], x[1]);
            t = temp[0];
            r = temp[1]
          })();
          jac = transform.jacobian_kruskal_to_schwarzschild(t, r);
          found_jac = (true);
        }
        if (chart == 1 && chart2 == 3) {
          (function() {
            var temp = [x[0], x[1]];
            t = temp[0];
            r = temp[1]
          })();
          jac = karl.array2d((2), (2));
          jac[0][0] = 1.0;
          jac[1][0] = 0.0;
          jac[0][1] = 0.0;
          jac[1][1] = (1.5) * Math.sqrt(r);
          found_jac = (true);
        }
        if (chart == 3 && chart2 == 1) {
          (function() {
            var temp = [x[0], x[1]];
            t = temp[0];
            u = temp[1]
          })();
          jac = karl.array2d((2), (2));
          jac[0][0] = 1.0;
          jac[1][0] = 0.0;
          jac[0][1] = 0.0;
          jac[1][1] = ((0.6666666666666666) * (Math.pow((u), (-0.3333333333333333))));
          found_jac = (true);
        }
        if (!found_jac) {
          throw io_util.strcat((["unrecognized charts, spacetime=", spacetime, ", chart=", chart, ", chart2=", chart2]));;
        }
        v2[0] = jac[0][0] * v[0] + jac[0][1] * v[1];
        v2[1] = jac[1][0] * v[0] + jac[1][1] * v[1];
        return (karl.clone_array1d(v2)); /* for python, divorce the new vector from entanglement with components of old */
      };
      transform.schwarzschild_to_kruskal = function(t, r) {
        var p, sign_b;

        /*
        Transforms Schwarzschild (t,r) coordinates to arcsinh-Kruskal (a,b). Assumes region I or II.
        */
        p = Math.sqrt(Math.abs(r - 1));
        if (r > 1.0) {
          sign_b = -1.0;
        } else {
          sign_b = 1.0;
        }
        return [transform.schwarzschild_to_kruskal_helper(p, r + t), sign_b * transform.schwarzschild_to_kruskal_helper(p, r - t)];
      };
      transform.schwarzschild_to_kruskal_helper = function(p, q) {

        /* compute y=asinh(pe^(q/2)) */
        return math_util.asinh_of_exp(Math.log(p) + q / 2);
      };
      transform.kruskal_to_schwarzschild = function(a, b) {
        var t, r, mu;

        /*
        Transforms arcsinh-Kruskal (a,b) coordinates to Schwarzschild (t,r).
        This is really just a wrapper for kruskal.aux().
        */
        (function() {
          var temp = kruskal.aux(a, b);
          t = temp[0];
          r = temp[1];
          mu = temp[2]
        })();
        return [t, r];
      };
      /*----------------------------------------------------------------------------------------------------- */
      transform.jacobian_schwarzschild_to_kruskal = function(t, r) {
        var jacobian, log_r_minus_1, xpi2, xmi2, s, q;

        /*
        Returns the matrix of partial derivatives of (a,b) with respect to (t,r), given t and r.
        The first index is 0 for a, 1 for b. The second index is 0 for t, 1 for r.
        The relevant equations are written up more legibly in the docs.
        For non-horizon points, assumes region I or II.
        For points on the horizon, the result contains some infinite matrix elements.
        */
        jacobian = karl.array2d((2), (2));
        log_r_minus_1 = Math.log(Math.abs(r - 1));
        xpi2 = math_util.safe_exp(-log_r_minus_1 - (r + t)); /* x_+^{-2} */
        xmi2 = math_util.safe_exp(-log_r_minus_1 - (r - t)); /* x_-^{-2} */
        if (r > 1.0) {
          s = 1.0; /* region I */
        } else {
          s = -1.0; /* region II */
        }
        jacobian[0][0] = 0.5 / Math.sqrt(1 + xpi2); /* da/dt */
        jacobian[1][0] = s * 0.5 / Math.sqrt(1 + xmi2); /* db/dt */
        if (r != 1.0) {
          /* not on the horizon */
          q = 1.0 / (1.0 - 1.0 / r);
          jacobian[0][1] = q * jacobian[0][0]; /* da/dr */
          jacobian[1][1] = -q * jacobian[1][0]; /* db/dr */
        } else {
          /* on the horizon */
          jacobian[0][1] = float("inf");
          jacobian[1][1] = float("inf");
        }
        return jacobian;
      };
      transform.jacobian_kruskal_to_schwarzschild = function(t, r) {
        var j1, a, d, b, c, det, j2;

        /*
        Returns the matrix of partial derivatives of (t,r) with respect to (a,b), given t and r.
        The first index is 0 for t, 1 for r. The second index is 0 for a, 1 for b.
        Slightly inefficient because we invert the 2x2 matrix,
        but not likely to affect performance because not called often.
        This can throw an error if called for a point on the horizon (r=1), where
        the Schwarzschild coordinates misbehave, and in general it's probably
        not going to be numerically accurate to use this near the horizon.
        */
        if (r == 1.0) {
          throw 'r=1 in transform.jacobian_kruskal_to_schwarzschild';;
        }
        j1 = transform.jacobian_schwarzschild_to_kruskal(t, r);
        a = j1[0][0];
        d = j1[1][1];
        b = j1[1][0];
        c = j1[0][1];
        det = a * d - b * c; /* should be nonzero because we checked above for r==1 */
        j2 = karl.array2d((2), (2)); /* allocate a new matrix, since a-d are just pointers into j1 */
        j2[1][1] = a / det;
        j2[0][0] = d / det;
        j2[1][0] = -b / det;
        j2[0][1] = -c / det;
        return j2;
      };
      /*----------------------------------------------------------------------------------------------------- */
      transform.schwarzschild_to_kruskal_small = function(t, r) {
        var t2, sc, cs, h, ks_t, ks_x, v, w;

        /*
        Compute Kruskal (a,b) coordinates from Schwarzschild (t,r). Overflows if t and r are not small.
        Returns a result in region I or II. For testing purposes only.
        */
        t2 = 0.5 * t;
        if (r > 1.0) {
          sc = Math.sinh(t2);
          cs = Math.cosh(t2);
        } else {
          sc = Math.cosh(t2);
          cs = Math.sinh(t2);
        }
        /* My formulation based on MTW p. 827: */
        h = Math.sqrt(Math.abs(r - 1.0)) * Math.exp(r / 2.0);
        ks_t = h * sc;
        ks_x = h * cs;
        v = ks_t + ks_x;
        w = ks_t - ks_x;
        return [Math.arcsinh(v), Math.arcsinh(w)];
      };
      transform.sch_is_in_future_light_cone = function(x, v) {
        var v_kruskal, va, vb, tol;

        /*
        Given vector v in Schwarzschild coordinates, determines whether it's in the future light cone.
        Returns [is_future,a_plus_b], where is_future is a boolean, and ab_sum is the value of v_a+v_b in
        KS coordinates. The boolean can be unreliable due to rounding if v is lightline; if you have a vector
        that is known lightlike and you just want to know whether it's future-oriented, it's better to test
        whether the sign of ab_sum is positive.
        */
        /* There is probably a more efficient way to do this, but this way seems bulletproof and manifestly correct. */
        v_kruskal = transform.transform_vector(v, x, 256, 1, {}, 2);
        va = v_kruskal[0];
        vb = v_kruskal[1];
        tol = 10 * (1.0e-16);
        return [(va >= -tol && vb >= -tol), va + vb];
      };
      transform.kruskal_to_time_zero = function(x, v, spacetime_or_chart, force) {
        var a, b, big, do_it, x2, v2, x3, v3, dt, did_it, x_s, v_s;

        /*
        Change kruskal coordinates to equivalent coordinates that correspond to Schwarzschild time t=0.
        If boolean input force is false, then nothing is done unless the point is the type of point near
        the horizon or photon sphere at large t for which this operation is likely to be helpful to precision.
        If the point is on the horizon, this operation is a no-op. As a convenience, this function can
        also be called when coords are not AKS, and then it's also a no-op.
        Returns [x2,v2,dt,did_it], where dt is the change in the Schwarzschild time that resulted from the
        transformation, and dit_it is true if anything actually happened.
        */
        if (spacetime_or_chart != (256 | 2)) {
          return [x, v, 0.0, (false)];
        }
        a = x[0];
        b = x[1];
        /* -- */
        big = 100.0;
        /* ...If this parameter is very big, like 10^6, then this code never gets executed and can't serve its */
        /*    purpose, while if I made it very small, like 1, then this code would get executed frequently, */
        /*    causing more rounding errors. */
        /*-- */
        do_it = (true);
        if (!force) {
          do_it = (Math.abs(a) > big * Math.abs(b)) || (Math.abs(b) > big * Math.abs(a));
        }
        if (!do_it) {
          return [x, v, 0.0, (false)];
        }
        if (b == 0.0 || a == 0.0) { /* can't define this operation for points on the horizon */
          return [x, v, 0.0, (false)];
        }
        /* If we're in region III or IV, then simply doing AKS->Schwarzschild->AKS would lose that information. */
        /* In that case, transform into I or II and recurse. */
        if (a < 0.0) {
          /* Invert in the origin in AKS. */
          x2 = (karl.clone_array1d(x));
          v2 = (karl.clone_array1d(v));
          x2[0] = -x2[0];
          x2[1] = -x2[1];
          v2[0] = -v2[0];
          v2[1] = -v2[1];
          (function() {
            var temp = transform.kruskal_to_time_zero(x2, v2, spacetime_or_chart, force);
            x3 = temp[0];
            v3 = temp[1];
            dt = temp[2];
            did_it = temp[3]
          })();
          x3[0] = -x3[0];
          x3[1] = -x3[1];
          v3[0] = -v3[0];
          v3[1] = -v3[1];
          return [x3, v3, dt, did_it];
        }
        /* Beyond this point, we're guaranteed to be in region I or II, and not at a horizon. */
        /* The following is a composition of three transformations: (1) AKS to Schwarzschild, (2) a time */
        /* translation, and (3) back to AKS. The jacobian for 2 is just the identity matrix, so the velocity */
        /* is not affected by it. */
        x_s = transform.transform_point(x, 256, 2, {}, 1);
        v_s = transform.transform_vector(v, x, 256, 2, {}, 1);
        dt = x_s[0];
        x_s[0] = 0.0;
        x2 = transform.transform_point(x_s, 256, 1, {}, 2);
        v2 = transform.transform_vector(v_s, x_s, 256, 1, {}, 2);
        return [x2, v2, dt, (true)];
      };
