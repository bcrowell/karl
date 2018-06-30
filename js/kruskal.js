      /*
               --- module kruskal ---
               This was translated from python. Do not edit directly.
            */
      if (typeof kruskal === 'undefined') {
        var kruskal = {};
      }

      /*
      Low-level routines to compute things for points in the Schwarzschild
      spacetime, in Kruskal-Szekeres coordinates, rescaled by an arcsinh,
      with the angular space embedded on the unit 2-sphere in a fictitious
      additional dimension.  The units are such that the Schwarzschild
      radius is 1, i.e., the mass is 1/2. 
      Transformations between arcsinh-Kruskal coordinates and Schwarzschild
      coordinates, and the associated jacobians, are in transform.pp, except
      that transformation from (a,b) to (t,r) are here, in kruskal.aux.
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
      karl.load("angular");
      karl.load("lambert_w_stuff");
      kruskal.metric_ks4 = function(p) {
        var g, t, r, mu, r2, sin_theta;

        /*
        For the Schwarzschild spacetime, compute the metric.
        Metric is in lower-index form, in +--- signature, with coordinates (a,b,theta,phi).
        */
        g = karl.array2d((4), (4));
        (function() {
          var temp = kruskal.aux(p[0], p[1]);
          t = temp[0];
          r = temp[1];
          mu = temp[2]
        })();
        g[0][1] = mu;
        g[1][0] = mu;
        r2 = r * r;
        sin_theta = Math.sin(p[2]);
        g[2][2] = -r2;
        g[3][3] = -r2 * sin_theta * sin_theta;
        return g;
      };
      kruskal.aux = function(a, b) {
        var e2a, e2b, f, u, r, mu, t;

        /*
        Compute Schwarzschild t and r, and metric element mu, for a point given in rescaled Kruskal coordinates (a,b).
        Return [t,r,mu].
        Although r is an analytic function of the coordinates throughout the entire
        maximal extension of the spacetime, numerical issues force us to do some case-splitting.
        The quantity mu is the coefficient in the line element, ds^2 = 2 mu da db -...,
        i.e., it's the matrix element g_ab of the metric. Returns t=None for points on a horizon.
        Returns None for all three variables if a and b lie beyond the singularity.
        There is also a C implementation of this function.
        */
        if (a < 0) {
          /* Flip to positive a in order to simplify some later computations. */
          return kruskal.aux(-a, -b);
        }
        /* From now on, we know a>=0. */
        e2a = math_util.safe_exp(-2 * a);
        e2b = math_util.safe_exp(-2 * Math.abs(b));
        /* From now on, we know we're in region I or II, a>0. */
        f = (1.0 - e2a) * (1.0 - e2b);
        u = a + Math.abs(b) + Math.log(f / 4.0) - 1.0;
        if (b < 0) {
          /* region I */
          r = 1.0 + lambert_w_stuff.lambert_w_of_exp(u);
        } else {
          /* region II */
          if (u > -1) {
            return [null, null, null]; /* beyond the singularity */
          }
          r = 1.0 + Math.lambert_w(-math_util.safe_exp(u));
        }
        /* Compute mu: */
        mu = (1.0 + e2a) * (1.0 + e2b) * (1 / (2 * Math.E * r)) * Math.exp(a + Math.abs(b) - (r - 1));
        /* Compute t: */
        if (a != 0 && b != 0 && (!(((r) == null) || (isNaN(r))))) {
          t = a - Math.abs(b) + Math.log((1 - e2a) / (1 - e2b));
        } else {
          t = null;
        }
        return [t, r, mu];
      };
      kruskal.metric = function(p) {
        var g, t, r, mu, r2;

        /*
        For the Schwarzschild spacetime, compute the metric.
        Metric is in lower-index form, in +--- signature, with coordinates (a,b,i,j,k).
        */
        g = karl.array2d((5), (5));
        (function() {
          var temp = kruskal.aux(p[0], p[1]);
          t = temp[0];
          r = temp[1];
          mu = temp[2]
        })();
        g[0][1] = mu;
        g[1][0] = mu;
        r2 = r * r;
        g[2][2] = -r2;
        g[3][3] = -r2;
        g[4][4] = -r2;
        return g;
      };
      kruskal.christoffel = function(p) {

        /*
        For the Schwarzschild spacetime, compute the Christoffel symbols, in coordinates (a,b,i,j,k).
        The order of indices is as in ctensor:
           symmetric on 1st 2 indices
           contravariant on final index
        See maxima/schwarzschild5.mac. In addition, we put in a fictitious centripetal term that
        keeps us constrained to the unit sphere i^2+j^2+k^2=1.
        */
        /*return christoffel_raw_maxima_output(p) */
        return kruskal.christoffel_massaged_maxima_output(p);
      };
      kruskal.christoffel_massaged_maxima_output = function(p) {
        var a, b, ch, t, r, mu, tanha, tanhb, q, q2, z;

        a = p[0];
        b = p[1];
        /* Based on christoffel_raw_maxima_output(), but simplified and massaged into a form that won't cause overflows. */
        ch = karl.array3d((5), (5), (5));
        /*------------------------------------------------------ */
        (function() {
          var temp = kruskal.aux(a, b);
          t = temp[0];
          r = temp[1];
          mu = temp[2]
        })();
        tanha = math_util.safe_tanh(a);
        tanhb = math_util.safe_tanh(b);
        q = (mu / 2) * (1.0 + 1.0 / r);
        q2 = -0.5 * mu / r;
        /*------------------------------------------------------ */
        ch[0][0][0] = tanha + q * tanhb; /* ^a_aa */
        ch[1][1][1] = tanhb + q * tanha;
        /*-- */
        z = q2 * tanhb;
        ch[0][2][2] = z; /* ^i_ai */
        ch[2][0][2] = z;
        ch[0][3][3] = z;
        ch[3][0][3] = z;
        ch[0][4][4] = z;
        ch[4][0][4] = z;
        /*-- */
        z = q2 * tanha;
        ch[1][2][2] = z; /* ^i_bi */
        ch[2][1][2] = z;
        ch[1][3][3] = z;
        ch[3][1][3] = z;
        ch[1][4][4] = z;
        ch[4][1][4] = z;
        /*-- */
        z = -0.5 * r * tanha; /* ^a_ii */
        ch[2][2][0] = z;
        ch[3][3][0] = z;
        ch[4][4][0] = z;
        /*-- */
        z = -0.5 * r * tanhb; /* ^b_ii */
        ch[2][2][1] = z;
        ch[3][3][1] = z;
        ch[4][4][1] = z;
        /*------------------------------------------------------ */
        angular.add_centripetal(ch, p);
        /*------------------------------------------------------ */
        return ch;
      };
      kruskal.christoffel_raw_maxima_output = function(p) {
        var a, b, ch;

        a = p[0];
        b = p[1];
        /* output of kruskal5.mac, plus centripetal terms */
        ch = karl.array3d((5), (5), (5));
        /*------------------------------------------------------ */
        ch[0][0][0] = ((Math.pow((((((Math.cosh((a))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))))) + (((2.0) * (Math.cosh((a))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.cosh((a))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))) * (((((2.0) * (Math.pow((Math.cosh((a))), (2.0))) * (Math.sinh((b))))) + (((Math.sinh((a))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))))) + (((((((Math.pow((Math.cosh((a))), (2.0))) * (Math.sinh((b))))) + (((2.0) * (Math.sinh((a))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.sinh((a))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))));
        /*   ... ^a _a a */
        ch[0][2][2] = ((-1.0) * (Math.cosh((a))) * (Math.sinh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^i _a i */
        ch[2][0][2] = ((-1.0) * (Math.cosh((a))) * (Math.sinh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^i _i a */
        ch[0][3][3] = ((-1.0) * (Math.cosh((a))) * (Math.sinh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^j _a j */
        ch[3][0][3] = ((-1.0) * (Math.cosh((a))) * (Math.sinh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^j _j a */
        ch[0][4][4] = ((-1.0) * (Math.cosh((a))) * (Math.sinh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^k _a k */
        ch[4][0][4] = ((-1.0) * (Math.cosh((a))) * (Math.sinh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^k _k a */
        ch[1][1][1] = ((Math.pow((((((Math.cosh((b))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))))) + (((2.0) * (Math.cosh((b))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.cosh((b))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))) * (((((2.0) * (Math.sinh((a))) * (Math.pow((Math.cosh((b))), (2.0))))) + (((Math.sinh((b))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))))) + (((((((Math.sinh((a))) * (Math.pow((Math.cosh((b))), (2.0))))) + (((2.0) * (Math.sinh((b))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.sinh((b))) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))));
        /*   ... ^b _b b */
        ch[1][2][2] = ((-1.0) * (Math.sinh((a))) * (Math.cosh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^i _b i */
        ch[2][1][2] = ((-1.0) * (Math.sinh((a))) * (Math.cosh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^i _i b */
        ch[1][3][3] = ((-1.0) * (Math.sinh((a))) * (Math.cosh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^j _b j */
        ch[3][1][3] = ((-1.0) * (Math.sinh((a))) * (Math.cosh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^j _j b */
        ch[1][4][4] = ((-1.0) * (Math.sinh((a))) * (Math.cosh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^k _b k */
        ch[4][1][4] = ((-1.0) * (Math.sinh((a))) * (Math.cosh((b))) * (Math.pow((((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) + (((2.0) * (Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))) + (((Math.pow((Math.E), (((1.0) + (Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))))))) * (Math.pow((Math.lambert_w((((-1.0) * (Math.pow((Math.E), (-1.0))) * (Math.sinh((a))) * (Math.sinh((b))))))), (2.0))))))), (-1.0))));
        /*   ... ^k _k b */
        ch[2][2][0] = -(Math.sinh(a) * Math.lambert_w(-Math.exp(-1) * Math.sinh(a) * Math.sinh(b)) + Math.sinh(a)) / (2 * Math.cosh(a));
        /*   ... ^a _i i */
        ch[2][2][1] = -(Math.sinh(b) * Math.lambert_w(-Math.exp(-1) * Math.sinh(a) * Math.sinh(b)) + Math.sinh(b)) / (2 * Math.cosh(b));
        /*   ... ^b _i i */
        ch[3][3][0] = -(Math.sinh(a) * Math.lambert_w(-Math.exp(-1) * Math.sinh(a) * Math.sinh(b)) + Math.sinh(a)) / (2 * Math.cosh(a));
        /*   ... ^a _j j */
        ch[3][3][1] = -(Math.sinh(b) * Math.lambert_w(-Math.exp(-1) * Math.sinh(a) * Math.sinh(b)) + Math.sinh(b)) / (2 * Math.cosh(b));
        /*   ... ^b _j j */
        ch[4][4][0] = -(Math.sinh(a) * Math.lambert_w(-Math.exp(-1) * Math.sinh(a) * Math.sinh(b)) + Math.sinh(a)) / (2 * Math.cosh(a));
        /*   ... ^a _k k */
        ch[4][4][1] = -(Math.sinh(b) * Math.lambert_w(-Math.exp(-1) * Math.sinh(a) * Math.sinh(b)) + Math.sinh(b)) / (2 * Math.cosh(b));
        /*   ... ^b _k k */
        /*------------------------------------------------------ */
        angular.add_centripetal(ch, p);
        /*------------------------------------------------------ */
        return ch;
      };
      kruskal.christoffel4 = function(p) {
        var ch, a, b, t, r, mu, ta, tb, ca, cb, sa, sb, q, big_b, eighth_r, sin_theta, sin2_theta, cos_theta;

        /*
        For the Schwarzschild spacetime, compute the Christoffel symbols, in coordinates (a,b,theta,phi).
        The order of indices is as in ctensor:
           symmetric on 1st 2 indices
           contravariant on final index
        */
        ch = karl.array3d((4), (4), (4));
        a = p[0];
        b = p[1];
        (function() {
          var temp = kruskal.aux(a, b);
          t = temp[0];
          r = temp[1];
          mu = temp[2]
        })();
        ta = math_util.safe_tanh(a);
        tb = math_util.safe_tanh(b);
        /* FIXME: In the following, we'll get overflows if a or b is large; massage the equations to avoid this problem. */
        ca = Math.cosh(a);
        cb = Math.cosh(b);
        sa = Math.sinh(a);
        sb = Math.sinh(b);
        q = (1.0 / r + 1.0 / (r * r));
        big_b = (4.0 / r) * Math.exp(-r);
        eighth_r = r / 8.0;
        sin_theta = Math.sin(theta);
        sin2_theta = Math.pow((sin_theta), (2.0));
        cos_theta = Math.cos(theta);
        ch[0][0][0] = ta + q * Math.exp(-r) * ca * sb; /* ^a_aa */
        ch[1][1][1] = tb + q * Math.exp(-r) * cb * sa; /* ^b_bb */
        ch[0][2][2] = -(big_b / (4 * r)) * sb / ca; /* ^theta_a theta */
        ch[2][0][2] = ch[0][2][2]; /* ^theta_theta a */
        ch[0][3][3] = ch[0][2][2]; /* ^phi_a phi */
        ch[3][0][3] = ch[0][2][2]; /* ^phi_phi a */
        ch[1][2][2] = -(big_b / (4 * r)) * sa / cb; /* ^theta_b theta */
        ch[1][3][3] = ch[1][2][2]; /* ^phi_b phi */
        ch[2][2][0] = -eighth_r * ta; /* ^a_theta theta */
        ch[2][2][1] = -eighth_r * tb; /* ^a_theta theta */
        ch[3][3][0] = -eighth_r * ta * sin2_theta; /* ^a_theta theta */
        ch[3][3][1] = -eighth_r * tb * sin2_theta; /* ^a_theta theta */
        ch[3][3][2] = -sin_theta * cos_theta; /* ^theta_phi phi */
        ch[2][3][3] = cos_theta / sin_theta; /* ^phi_theta phi */
        ch[3][2][3] = ch[2][3][3]; /* ^phi_phi theta */
        return ch;
      };
