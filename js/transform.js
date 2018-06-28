      /*
               --- module transform ---
               This was translated from python. Do not edit directly.
            */
      if (typeof transform === 'undefined') {
        var transform = {};
      }

      /*
      Low-level routines to compute transformations between Schwarzschild and
      arcsinh-Kruskal coordinates, and the jacobians of those transformations.
      Documentation for the math is in the file doc.tex, which can be
      compiled to pdf format by doing a "make doc." (Comments in the code do
      not document the math or the definitions of the variables.) 
      */
      /* ... note that (NaN)==(NaN) is false, so use IS_(NaN) */
      karl.load("lib/array");;
      /* ... works in rhino and d8 */
      /* ... https://stackoverflow.com/q/26738943/1142217 */
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
      karl.load("math_util");
      karl.load("kruskal");
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
