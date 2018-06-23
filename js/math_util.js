      /*
               --- module math_util ---
               This was translated from python. Do not edit directly.
            */
      if (typeof math_util === 'undefined') {
          var math_util = {};
      }

      ;

      ;;


      /* ... relative precision for arithmetic */
      /* ... ln of (1.0e-16) */
      /* ... ln of greatest number we can store in floating point; IEEE-754 floating point can store 2^128-1 */
      math_util.safe_exp = function(x) {

          /*
          Return exp(x), or 0 if x is so large and negative that exp(x) is negligible compared to unity,
          and might have caused an underflow.
          */
          if (x < (-36.8413614879047)) {
              return 0.0;
          }
          return Math.exp(x);
      };
      math_util.safe_tanh = function(x) {
          var z;

          /*
          Return tanh(x), without overflows, underflows, or unnecessarily inefficient evaluation of negligible stuff.
          */
          if (x < 0.0) {
              return -math_util.safe_tanh(-x);
          }
          /* From here on, we know x>=0. */
          z = math_util.safe_exp(-2 * x);
          return (1.0 - z) / (1.0 + z);
      };
      math_util.asinh_of_exp = function(u) {
          var y, g, m, coeffs, term;

          /*
          Compute asinh(e^u), using asymptotics if u is large.
          Series expansion: https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Series_expansions
          */
          if (u < 0.5 * (88.722839111673)) {
              /* ... IEEE-754 floating point can store exp(88); asinh probably uses the square of the input internally */
              return Math.arcsinh(Math.exp(u));
          }
          y = Math.log(2.0) + u;
          if (u > (-36.8413614879047)) {
              return y; /* correction term will be negligible compared to unity */
          }
          g = Math.exp(-2.0 * u); /* 1/x^2, where x is the input to the arcsinh */
          m = 1.0;
          coeffs = [1.0 / 4.0, -3.0 / 32.0, 15.0 / 288.0, -105.0 / 3072.0];
          for (var n = 1; n < 6; n++) {
              m = m * g;
              term = coeffs[n - 1] * m;
              /* print("calculating correction, n=",n,", term=",term) */
              if (term < (1.0e-16)) {
                  return y;
              }
              y = y + term;
          }
          return y; /* If we get here, the approx. may not be very good. */
      };
