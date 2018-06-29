// (c) 2018 Benjamin Crowell, GPL v3

(function() {

  var elliptical_orbit_period = function(r, a, direction, n) {
    var spacetime, chart, v_phi, x, circular_period, v, q, r_max, period, delta_lambda, ndebug, opt, err, final_x, final_v, final_lambda, info, eps;

    /*
    Period of an elliptical orbit, Schwarzschild coordinates.
    Start at perihelion, r. Make the initial velocity greater than the circular-orbit value by the factor a.
    direction = angle about the x axis for the initial motion, defines plane of orbit
    Runge-Kutta with n steps.
    */
    spacetime = SP_SCH; // Schwarzschild spacetime
    chart = CH_SCH; // sch5 coordinates
    v_phi = 1 / Math.sqrt(2.0 * r * r * r); /* exact condition for circular orbit in Sch., if v_t=1. */
    x = [0.0, r, 1.0, 0.0, 0.0];
    circular_period = 2.0 * Math.PI / v_phi;
    /* Increase velocity at perihelion to make orbit elliptical: */
    v_phi = v_phi * a;
    v = [1.0, 0.0, 0.0, v_phi * Math.cos(direction), v_phi * Math.sin(direction)];
    v = angular.make_tangent(x, v);
    v = vector.normalize(spacetime, chart, x, v);
    /* Compute newtonian r_max: */
    q = a * a / 2.0 - 1
    r_max = r * ((-1.0 - Math.sqrt(1.0 + 2.0 * a * a * q)) / (2 * q))
    period = circular_period * Math.pow((r + r_max) / (r + r), 1.5) // Kepler's law of periods
    delta_lambda = period;
    /*-- */
    ndebug = 0;
    var nchunks = 10; // n should be divisible by this
    var delta_lambda_per_chunk = delta_lambda/nchunks;
    opt = {
      'dlambda': delta_lambda / n,
      'ndebug': ndebug,
      'norm_final': false
    };
    lam = 0.0;
    var state = {"chunk":0,"data":[lam,x,v],"claimed":0};
      // state.claimed=0 means it's ready for someone to claim it, 1 means someone's working on it
    var finish_up = function(err,message) {
      clearInterval(interval_id);
      print("final x=", io_util.vector_to_str_n_decimals(state.data[1], 16));
      if (err) {print(message)}
      print("...done");
    };
    var output_helper = function(i,lam,x,v) {
      print("step=",i," lambda/T=",io_util.fl(lam/period), 
                      " x=(",io_util.vector_to_str_n_decimals(x,2), 
                      ") v=(",io_util.vector_to_str_n_decimals(v,2),")");
    };
    var run_some_steps = function() {
      if (state.claimed==1) {return;}
      state.claimed=1;
      var chunk = state.chunk;
      var lam = state.data[0];
      var x   = state.data[1];
      var v   = state.data[2];
      opt.lambda0 = lam;
      opt.lambda_max = lam+delta_lambda_per_chunk;
      var temp = runge_kutta.trajectory_simple(spacetime, chart, x, v, opt); // returns [0,x,v,lam,{}]
      err = temp[0];
      x = temp[1];
      v = temp[2];
      lam = temp[4];
      info = temp[5];
      chunk += 1;
      output_helper(n*chunk/nchunks,lam,x,v);
      state.chunk = chunk;
      state.data[0] = lam;
      state.data[1] = x;
      state.data[2] = v;
      if (chunk>=nchunks || err) {
        finish_up(err,info["message"]);
      }
      state.claimed=0;
    }
    print("T=",io_util.fl(period));
    output_helper(0,lam,x,v);
    var interval_id = setInterval(run_some_steps,10); // time in milliseconds
  };

  var r, a, direction, n;
  r = 1.0e8;
  a = 1.1;
  direction = 0.0;
  n = 100000;
  print("starting...");
  elliptical_orbit_period(r, a, direction, n);

})();
