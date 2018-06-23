      /*
               --- module io_util ---
               This was translated from python. Do not edit directly.
            */
      if (typeof io_util === 'undefined') {
          var io_util = {};
      }

      ;
      /* l is an array of objects which may not be strings */
      io_util.strcat = function(l) {
          return l.join("");
      };
      io_util.print_no_newline = function(s) {

          PRINT(s);
      };
      io_util.vector_to_str = function(v) {

          return io_util.vector_to_str_n_decimals(v, 3);
      };
      io_util.vector_to_str_n_decimals = function(v, n) {
          var f = [];
          for (var i = 0; i < v.length; i++) {
              f.push(v[i].toExponential(n).toString());
          }
          return f.join(',');
      }
      io_util.fl = function(x) {
          return x.toExponential(3).toString();
      }
