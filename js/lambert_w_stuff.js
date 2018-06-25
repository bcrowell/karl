      /*
               --- module lambert_w_stuff ---
               This was translated from python. Do not edit directly.
            */
      if (typeof lambert_w_stuff === 'undefined') {
        var lambert_w_stuff = {};
      }


      ;;;



      /* ... note that (NaN)==(NaN) is false, so use IS_(NaN) */
      karl.load("lib/array");;
      /* ... works in rhino and d8 */
      /* ... https://stackoverflow.com/q/26738943/1142217 */
      /*           ... see notes above about usage with array literals */
      /*           ... in JS, numbers are primitives, not objects, so no need clone them */
      /*################################################################## */
      /*################################################################## */
      lambert_w_stuff.lambert_w_of_exp = function(u) {
        var l2, nterms, w, niter, z, q, err;

        /*
        Return W(e^u), where W is the Lambert W function W0.
        The purpose of this function is to handle cases where u is too big
        to allow x=e^u to be stored in floating point. The method used is the
        iterative scheme described in Veberic, https://arxiv.org/abs/1003.1628 , 
        sec. 2.3. The output W of the function satisfies u=ln W+W to within a
        relative error of 10^-16.
        */
        if (u < 100) {
          return Math.lambert_w(Math.exp(u));
        }
        /* Find an initial guess, https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions . */
        l2 = Math.log(u); /* ln(ln(x)) */
        nterms = Math.floor(88 / l2); /* truncate series to avoid underflow; IEEE-754 can represent 2^128; here 88=ln(2^128) */
        w = u - l2;
        if (nterms >= 2) {
          w += l2 / u;
          if (nterms >= 3) {
            w += 0.5 * l2 * (-2 + l2) / (u * u);
            if (nterms >= 4) {
              w += l2 * (6 - 9 * l2 + 2 * l2 * l2) / (6 * u * u * u);
            }
          }
        }
        /* Iteration: */
        /* Predetermine how many iterations we need, based on empirical tests. */
        if (nterms >= 5) {
          niter = 2;
        } else {
          niter = 1;
        }
        for (var i = 0; i < niter; i++) {
          z = u - Math.log(w) - w;
          q = 2 * (1 + w) * (1 + w + 2.0 * z / 3.0);
          err = (z / (1 + w)) * ((q - z) / (q - 2 * z));
          w = w * (1 + err);
        }
        return w;
      };
