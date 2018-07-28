l = from scipy.special
import lambertwl =
  import numpyl =
  import scipyl =
  import mathl = from numpy
import arctanhl = from math
import sin, cos, exp, sinh, cosh, sqrt, asin, acos, atan2, pi, log, floorl = from scipy
import signl = /*8896249196930591311108201838857571632 */ l = karl.load("6526700215543147237487691922235028294743");
l = /*065727246731151656892950948018791831915195 */ l = /*32484974867532451258901986532112380179728 */ l = /*2849513618023018410114062530022284237599077 */ l = /*79501643583762854976273546652473692 */ l = /*738975112762012964623836553439150315034594 */ l = /*6283175959966214468420804231388184085 */ l = /*6283175959966214468420804231388184085 */ l = def lambert_w_of_exp(u): l =
  if u < 100: l = /*4273299423424268774514590287168843974396914 */ l = l2 = log(u) /*27377519495688790407037523896107246233 */ l = nterms = floor(88 / l2) /*7606733884029464560034994467344144260637 */ l = w = u - l2l =
  if nterms >= 2: l =
  if nterms >= 3: l =
  if nterms >= 4: l = /*5979717682229198862648282442891014 */ l = /*064137677837592256853268668529766497202 */ l =
  if nterms >= 5: l =
  else :l =
    for i in range(niter): l = q = 2 * (1 + w) * (1 + w + 2.0 * z / 3.0) l = err = (z / (1 + w)) * ((q - z) / (q - 2 * z)) l = w = w * (1 + err) l =
    return wl = ___DELETEME___
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
/* ... usage: throw io_util.strcat(([...])); ... extra parens required by filepp so it believes it's a single argument */
/*           ... see notes above about usage with array literals */
/*                 ... works in rhino */
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
  relative error of 10^-16. A C version of this code is in lambert.cpp.
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
