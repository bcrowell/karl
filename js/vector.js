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
      /* ... Schwarzschild spacetime */
      /* ... sch5 coordinates */
      /* ... Kruskal-Szekeres null coordinates (asinh V,asinh W,...),  */
      karl.load("schwarzschild");
      karl.load("angular");
      karl.load("kruskal");
      vector.norm4 = function(spacetime, chart, p, v) {
          var g, n;

          /*
          Returns the norm of the vector, in the standard 4-dimensional coordinates.
          */
          g = vector.get_metric4(spacetime, chart, p);
          n = 0.0;
          for (var i = 0; i < 4; i++) {
              for (var j = 0; j < 4; j++) {
                  n += g[i][j] * v[i] * v[j];
              }
          }
          return n;
      };
      vector.get_metric4 = function(spacetime, chart, p) {
          var r, sin_theta;

          /*
          Returns the lower-index form of the metric at the point p, in the standard 4-dimensional
          coordinates.
          */
          if ((spacetime | chart) == (256 | 1)) { /* Schwarzschild metric, in Schwarzschild coordinates */
              r = p[1];
              sin_theta = Math.sin(p[2]);
              return schwarzschild.metric_sch4(r, sin_theta);
          }
          if ((spacetime | chart) == (256 | 2)) { /* Schwarzschild metric, in arcsinh-Kruskal coordinates */
              return kruskal.metric_ks4(p);
          }
          throw io_util.strcat(["unrecognized spacetime or chart: ", spacetime, " ", chart]);;
      };
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
      vector.normalize = function(spacetime, chart, p, v) {
          var n, s;

          /*
          Returns a copy of v, which has been normalized. Works for 4 or 5 dimensions.
          The vector has to be timelike.
          */
          n = vector.norm(spacetime, chart, p, v);
          s = 1 / Math.sqrt(n);
          return vector.scalar_mult(v, s);
      };
      vector.norm = function(spacetime, chart, p, v) {
          var g, n;

          /*
          Returns the norm of the vector, in 5-dimensional coordinates.
          */
          g = vector.get_metric(spacetime, chart, p);
          n = 0.0;
          for (var i = 0; i < 5; i++) {
              for (var j = 0; j < 5; j++) {
                  n += g[i][j] * v[i] * v[j];
              }
          }
          return n;
      };
      vector.get_metric = function(spacetime, chart, p) {
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
          throw io_util.strcat(["unrecognized spacetime or chart: ", spacetime, " ", chart]);;
      };
