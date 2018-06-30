      /*
               --- module angular ---
               This was translated from python. Do not edit directly.
            */
      if (typeof angular === 'undefined') {
        var angular = {};
      }

      /*
      Helper routines to compute things for the representation of the
      unit sphere as a sphere embedded in three-dimensional cartesian
      space (i,j,k).
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
      angular.renormalize = function(x) {
        var n, xi;

        /*
        Given a point represented by coordinates of the form (...,...,i,j,k),
        return a copy with i, j, and k renormalized to lie on the unit sphere.
        */
        n = (karl.clone_array1d(x));
        xi = Math.sqrt(n[2] * n[2] + n[3] * n[3] + n[4] * n[4]);
        n[2] = n[2] / xi;
        n[3] = n[3] / xi;
        n[4] = n[4] / xi;
        return n;
      };
      angular.make_tangent = function(x, v0) {
        var v, dot;

        /*
        Given a point x and a velocity vector v in its tangent space, return a copy of
        v in which the component perpendicular to the i-j-k sphere has been removed.
        For accurate results, x should already be accurately normalized. Since only the
        angular parts of x and v0 are used, it doesn't matter whether they are expressed
        in Schwarzschild or Kruskal, as long as they're 5-dimensional.
        */
        v = (karl.clone_array1d(v0));
        dot = x[2] * v[2] + x[3] * v[3] + x[4] * v[4];
        v[2] -= dot * x[2];
        v[3] -= dot * x[3];
        v[4] -= dot * x[4];
        return v;
      };
      angular.theta_phi_to_ijk = function(theta, phi) {
        var sin_theta;

        sin_theta = Math.sin(theta);
        return [sin_theta * Math.cos(phi), sin_theta * Math.sin(phi), Math.cos(theta)];
      };
      angular.add_centripetal = function(ch, p) {
        var i, j, k, xi2, z;

        /*
        Add centripetal terms to the Christoffel symbols for a 5-dimensional coordinate system.
        Modifies ch in place by adding the centripetal parts.
        */
        i = p[2];
        j = p[3];
        k = p[4];
        xi2 = i * i + j * j + k * k;; /* should normally be very close to 1 */
        for (var m = 2; m < 5; m++) {
          z = p[m];
          for (var n = 2; n < 5; n++) {
            ch[n][n][m] += z / xi2;
          }
        }
      };
