extern "C" {
  double veberic_lambert_w(double);
  double lambert_w_of_exp(double);
}
#include <math.h>

// Cut-down version of Veberic's LambertW.h:
namespace utl {
  template<int branch>
  double LambertW(double);
}

// C wrapper for Veberic's C++ function/template setup, which I don't understand.
double veberic_lambert_w(double x) {
  return utl::LambertW<0>(x);
}

// C implementation of the function of the same name in lambert_w_stuff. See
// comments there.
double lambert_w_of_exp(double u) {
  double l2,w,z,q,err;
  int i,nterms,niter;
  if (u<100.0) {
    return veberic_lambert_w(exp(u));
  }
  // Find an initial guess, https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions .
  l2 = log(u); // =ln(ln(x))
  nterms = floor(88.0/l2);
  // ... truncate series to avoid underflow; IEEE-754 can represent 2^128; here 88=ln(2^128)
  w = u-l2;
  if (nterms>=2) {
    w += l2/u;
    if (nterms>=3) {
      w += 0.5*l2*(-2.0+l2)/(u*u);
      if (nterms>=4) {
        w += l2*(6-9*l2+2*l2*l2)/(6*u*u*u);
      }
    }
  }
  // Iteration:
  // Predetermine how many iterations we need, based on empirical tests.
  if (nterms>=5) {niter=2;} else {niter=1;}
  for (i=0; i<niter; i++) {
    z = u-log(w)-w;
    q = 2*(1+w)*(1+w+2.0*z/3.0);
    err = (z/(1+w))*((q-z)/(q-2*z));
    w = w*(1+err);
  }
  return w;
}
