      /*
               --- module vector ---
               This was translated from python. Do not edit directly.
            */
      if (typeof vector === 'undefined') {
        var vector = {};
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
      karl.load("schwarzschild");
      karl.load("angular");
      karl.load("kruskal");
      karl.load("io_util");
      vector.scalar_mult = function(v0, s) {
        var v, i;

        /*
        Returns a copy of v, which has been multiplied by the scalar s. Works for 4 or 5 dimensions.
        */
        v = (karl.clone_array1d(v0));
        for (var i = 0; i < v.length; i++) {
          v[i] = v[i] * s;
        }
        return v;
      };
      vector.negate = function(v) {

        return vector.scalar_mult(v, -1.0);
      };
      vector.add = function(u, v) {
        var w, i;

        /*
        Returns the sum of vectors u and v.
        */
        w = (karl.clone_array1d(u));
        for (var i = 0; i < u.length; i++) {
          w[i] = w[i] + v[i];
        }
        return w;
      };
      vector.sub = function(u, v) {
        var w, i;

        /*
        Returns u-v.
        */
        w = (karl.clone_array1d(u));
        for (var i = 0; i < u.length; i++) {
          w[i] = w[i] - v[i];
        }
        return w;
      };
      vector.normalize = function(spacetime, chart, pars, p, v) {
        var n, s;

        /*
        Returns a copy of v, which has been normalized. Works for 4 or 5 dimensions.
        The vector has to be timelike. To normalize a spacelike vector, see normalize_spacelike().
        */
        n = vector.norm(spacetime, chart, pars, p, v);
        s = 1 / Math.sqrt(n);
        return vector.scalar_mult(v, s);
      };
      vector.normalize_spacelike = function(spacetime, chart, pars, p, v) {
        var n, s;

        /*
        Returns a copy of v, which has been normalized so that |v|^2=-1. Works for 4 or 5 dimensions.
        The vector has to be spacelike.
        */
        n = vector.norm(spacetime, chart, pars, p, v);
        s = 1 / Math.sqrt(-n);
        return vector.scalar_mult(v, s);
      };
      vector.norm = function(spacetime, chart, pars, p, v) {
        var g, n;

        /*
        Returns the squared norm of the vector, which should be supplied in 5-dimensional coordinates.
        */
        g = vector.get_metric(spacetime, chart, pars, p);
        n = 0.0;
        for (var i = 0; i < 5; i++) {
          for (var j = 0; j < 5; j++) {
            n += g[i][j] * v[i] * v[j];
          }
        }
        return n;
      };
      vector.inner_product = function(spacetime, chart, pars, p, u, v) {
        var g, result;

        /*
        Returns the norm of the vector, in 5-dimensional coordinates.
        */
        g = vector.get_metric(spacetime, chart, pars, p);
        result = 0.0;
        for (var i = 0; i < 5; i++) {
          for (var j = 0; j < 5; j++) {
            result += g[i][j] * u[i] * v[j];
          }
        }
        return result;
      };
      vector.proj = function(spacetime, chart, pars, p, u, v) {
        var s;

        /*
        Returns P_u(v)=v-[(v.u)/(u.u)]u, the part of v that is perpendicular to u.
        The vector u should be future timelike.
        */
        s = vector.inner_product(spacetime, chart, pars, p, u, v) / vector.norm(spacetime, chart, pars, p, u);
        return vector.add(v, vector.scalar_mult(u, -s));
      };
      vector.proj_spacelike = function(spacetime, chart, pars, p, u, v) {
        var s;

        /*
        Returns v+[(v.u)/(u.u)]u, the part of v that is perpendicular to u.
        Both u and v should be spacelike.
        */
        s = vector.inner_product(spacetime, chart, pars, p, u, v) / vector.norm(spacetime, chart, pars, p, u);
        return vector.add(v, vector.scalar_mult(u, s));
      };
      vector.get_metric = function(spacetime, chart, pars, p) {
        var r;

        /*
        Returns the lower-index form of the metric at the point p, in 5-dimensional
        coordinates.
        */
        if ((spacetime | chart) == (256 | 1)) { /* Schwarzschild metric, in Schwarzschild coordinates */
          r = p[1];
          return schwarzschild.metric(r);
        }
        if ((spacetime | chart) == (256 | 2)) { /* Schwarzschild metric, in arcsinh-Kruskal coordinates */
          return kruskal.metric(p);
        }
        throw io_util.strcat((["unrecognized spacetime or chart: ", spacetime, " ", chart]));;
      };
