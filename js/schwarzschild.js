      /*
               --- module schwarzschild ---
               This was translated from python. Do not edit directly.
            */
      if (typeof schwarzschild === 'undefined') {
        var schwarzschild = {};
      }

      /*
      Low-level routines to compute things for points in the Schwarzschild
      spacetime, in Schwarzschild coordinates, with the angular space
      embedded on the unit 2-sphere in a fictitious additional dimension.
      The units are such that the Schwarzschild radius is 1, i.e., the mass
      is 1/2.  This file does not contain any of the following
      functionality: Kruskal coordinates, coordinate transformations, operations 
      on vectors.
      Some routines in this code compute things relating to the the normal
      4-dimensional Schwarzschild coordinates, and these have named with
      sch4 in them.
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
      if (!(typeof window !== 'undefined') && (typeof Math.karl === 'undefined')) {
        /* load() works in rhino,  !  sure about other engines */
        load("lib/math.js");
        if (typeof one_over_E === 'undefined') {
          load("lib/lambertw.js")
        }
        load("lambert_w_stuff.js");
        Math.karl = 1;
      }
      schwarzschild.metric_sch4 = function(r, sin_theta) {
        var g, a, r2;

        /*
        For the Schwarzschild spacetime, compute the metric, in Schwarzschild coordinates.
        The metric is in lower-index form, in +--- signature, with coordinates (t,r,theta,phi).
        The mass is assumed to be 1/2, so that r is in units of the Schwarzschild radius.
        Angles in radians.
        Metric taken from https://en.wikipedia.org/wiki/Schwarzschild_metric .
        */
        g = karl.array2d((4), (4));
        a = 1.0 - 1.0 / r;
        r2 = r * r;
        g[0][0] = a;
        g[1][1] = -1.0 / a;
        g[2][2] = -r2;
        g[3][3] = -r2 * sin_theta * sin_theta;
        return g;
      };
      schwarzschild.metric = function(r) {
        var g, a, r2;

        /*
        For the Schwarzschild spacetime, compute the metric, in 5-dimensional Schwarzschild coordinates.
        The metric is in lower-index form, in +--- signature, with coordinates (t,r,i,j,k).
        The mass is assumed to be 1/2, so that r is in units of the Schwarzschild radius.
        Angles in radians.
        */
        g = karl.array2d((5), (5));
        a = 1.0 - 1.0 / r;
        r2 = r * r;
        g[0][0] = a;
        g[1][1] = -1.0 / a;
        g[2][2] = -r2;
        g[3][3] = -r2;
        g[4][4] = -r2;
        return g;
      };
      schwarzschild.christoffel = function(p) {
        var ch, r, r2, c, z, n, i, j, k, xi2, m;

        /*
        For the Schwarzschild spacetime, compute the Christoffel symbols, in Schwarzschild coordinates (t,r,i,j,k).
        The order of indices is as in ctensor:
           symmetric on 1st 2 indices
           contravariant on final index
        See maxima/schwarzschild5.mac. In addition, we put in a fictitious centripetal term that
        keeps us constrained to the unit sphere i^2+j^2+k^2=1.
        */
        ch = karl.array3d((5), (5), (5));
        /* The following symbols, involving only t and r, are the same as the standard ones in sch4: */
        r = p[1];
        r2 = r * r;
        c = r - 1.0; /* = r-2GM in Carroll's notation */
        ch[0][0][1] = 0.5 * c / (r2 * r);
        z = 0.5 / (r * c);
        ch[1][1][1] = -z;
        ch[0][1][0] = z;
        ch[1][0][0] = z;
        /* These terms are analogous to Gamma^r_theta_theta=-c and Gamma^theta_r_theta=1/r in sch4: */
        for (var n = 2; n < 5; n++) {
          ch[n][n][1] = -c;
          ch[n][1][n] = 1 / r;
          ch[1][n][n] = 1 / r;
        }
        /* Fictitious centripetal terms: */
        i = p[2];
        j = p[3];
        k = p[4];
        xi2 = i * i + j * j + k * k;; /* should normally be very close to 1 */
        for (var m = 2; m < 5; m++) {
          z = p[m];
          for (var n = 2; n < 5; n++) {
            ch[n][n][m] = z / xi2;
          }
        }
        return ch;
      };
      schwarzschild.christoffel_sch4 = function(t, r, sin_theta, cos_theta) {
        var ch, r2, c, z, s2, cot;

        /*
        For the Schwarzschild spacetime, compute the Christoffel symbols, in Schwarzschild coordinates (t,r,theta,phi).
        The order of indices is as in ctensor:
           symmetric on 1st 2 indices
           contravariant on final index
        This will give a division by zero exception if used at the coordinate singularities at theta=0 and pi
        and r=1.
        The equations are from http://ned.ipac.caltech.edu/level5/March01/Carroll3/Carroll7.html ,
        equation 7.33.
        */
        ch = karl.array3d((4), (4), (4));
        r2 = r * r;
        c = r - 1.0; /* = r-2GM in Carroll's notation */
        ch[0][0][1] = 0.5 * c / (r2 * r);
        z = 0.5 / (r * c);
        ch[1][1][1] = -z;
        ch[0][1][0] = z;
        ch[1][0][0] = z;
        z = 1.0 / r;
        ch[1][2][2] = z;
        ch[2][1][2] = z;
        ch[2][2][1] = -c;
        ch[1][3][3] = z;
        ch[3][1][3] = z;
        s2 = sin_theta * sin_theta;
        ch[3][3][1] = -c * s2;
        ch[3][3][2] = -sin_theta * cos_theta;
        cot = cos_theta / sin_theta;
        ch[2][3][3] = cot;
        ch[3][2][3] = cot;
        return ch;
      };
