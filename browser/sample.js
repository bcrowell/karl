// (c) 2018 Benjamin Crowell, GPL v3

(function() {

  function die(message) {
    throw new Error(message);
    stop_animation(); // Why doesn't throwing the error already terminate the program?
  }
  function display_message(m) {
    document.getElementById("message").innerHTML = m;
  }


  var elliptical_orbit_period = function(r, a, direction, n) {
  var spacetime, chart, v_phi, x, circular_period, v, q, r_max, period, lambda_max, ndebug, opt, err, final_x, final_v, final_lambda, info, eps;

  /*
  Period of an elliptical orbit, Schwarzschild coordinates.
  Start at perihelion, r. Make the initial velocity greater than the circular-orbit value by the factor a.
  Test against the Keplerian period. There is no point in testing with large n, because the errors
  become dominated by the Keplerian approximation.
  direction = angle about the x axis for the initial motion, defines plane of orbit
  Runge-Kutta with n steps.
  This is intended to be used with very large r, so that the Keplerian approximation is good.
  */
  spacetime = 256;
  chart = 1;
  v_phi = 1 / Math.sqrt(2.0 * r * r * r); /* exact condition for circular orbit in Sch., if v_t=1. */
  x = [0.0, r, 1.0, 0.0, 0.0];
  circular_period = 2.0 * Math.PI / v_phi;
  /* Increase velocity at perihelion to make orbit elliptical: */
  v_phi = v_phi * a;
  v = [1.0, 0.0, 0.0, v_phi * Math.cos(direction), v_phi * Math.sin(direction)];
  v = angular.make_tangent(x, v);
  v = vector.normalize(spacetime, chart, x, v);
  /* Compute newtonian r_max: */
  q = ((-1.0) + (((0.5) * (Math.pow((a), (2.0))))));
  r_max = ((((1.0) / (2.0))) * (Math.pow((q), (-1.0))) * (((-1.0) + (((-1.0) * (Math.pow((((1.0) + (((2.0) * (Math.pow((a), (2.0))) * (q))))), (((1.0) / (2.0))))))))) * (r));
  period = ((0.35355339059327384) * (circular_period) * (Math.pow((((Math.pow((r), (-1.0))) * (((r) + (r_max))))), (1.5)))); /* Kepler's law of periods */
  lambda_max = period;
  /*-- */
  ndebug = 0;
  if (verbosity >= 3) {
      ndebug = n / 10;
  }
  opt = {
      'lambda_max': lambda_max,
      'dlambda': lambda_max / n,
      'ndebug': ndebug,
      'norm_final': (false)
  };
  (function() {
      var temp = runge_kutta.geodesic_simple(spacetime, chart, x, v, opt);
      err = temp[0];
      final_x = temp[1];
      final_v = temp[2];
      final_lambda = temp[3];
      info = temp[4]
  })();
  if (verbosity >= 2) {
      print("final x=", io_util.vector_to_str_n_decimals(final_x, 16));
  }
  };

  var r,a,direction,n;
      r = 1.0e8;
      a = 1.1;
      direction = 0.0;
      verbosity=3; // global variable
      n = 100;
  print("starting...");
  elliptical_orbit_period(r,a,direction,n);
  print("...done");

})();
