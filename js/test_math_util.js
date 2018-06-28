      /*
               --- main program test_math_util ---
               This was translated from python. Do not edit directly.
            */
      var karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
        load("lib/loader.js");
      }
      if (typeof test_math_util === 'undefined') {
        var test_math_util = {};
      }
      var assert_rel_equal, assert_equal, assert_rel_equal_eps, assert_equal_eps, asinh_of_exp_tests, u, v, y, yy;

      /*!/usr/bin/python3 */
      /* ... note that (NaN)==(NaN) is false, so use IS_(NaN) */
      karl.load("lib/array");;
      /* ... works in rhino and d8 */
      /* ... https://stackoverflow.com/q/26738943/1142217 */
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
      /* This should only be included (and executed) once, by the main program. */
      load("test.js");;
      assert_rel_equal = test.assert_rel_equal;;
      assert_equal = test.assert_equal;;
      assert_rel_equal_eps = test.assert_rel_equal_eps;;
      assert_equal_eps = test.assert_equal_eps;;
      /* ... relative precision for arithmetic */
      /* ... ln of (1.0e-16) */
      /* ... ln of greatest number we can store in floating point; IEEE-754 floating point can store 2^128-1 */
      karl.load("math_util");
      /* ------------Test safe_exp(). -------------------- */
      test.assert_rel_equal_eps(2.71828182845904523536028747135266249775724709369995, math_util.safe_exp(1.0), (1.0e-16) * 5);
      test.assert_rel_equal_eps(0.0, math_util.safe_exp(-1.0e50), (1.0e-16) * 5);
      /* ------------Test safe_tanh(). -------------------- */
      test.assert_rel_equal_eps(0.761594155955765, math_util.safe_tanh(1.0), (1.0e-16) * 20);
      test.assert_rel_equal_eps(1.0, math_util.safe_tanh(1.0e50), (1.0e-16) * 5);
      /* ------------Test asinh_of_exp_tests(). -------------------- */
      /* The following was calculated by misc/generate_asinh_of_exp_tests.py : */
      asinh_of_exp_tests = [
        [1, 0.032235378292369784527865848946712188199057140312684],
        [6, 0.0000015360495491554979773791726417095636553617222490443],
        [11, 0.000000000069736702314428308718040371746639796885029839817346],
        [16, 0.0000000000000031660413872735288950531271456171885279637586456956],
        [21, 1.4373805660733899513561857863845093903996543240323e-19],
        [26, 6.5256976741692620117566744046126154463510038869117e-24],
        [31, 2.9626621605849525157126680815586589786020927176369e-28],
        [36, 1.3450465400052846075367888954673816013582613197243e-32],
        [41, 6.1065018443509058234726435961025347510901368590978e-37],
        [46, 2.7723475477407621994868503415523225432864414014397e-41],
        [51, 1.2586730574825222807565552060615451795219859983619e-45],
        [56, 3.6082323586244641222897242306545005364005157881418e-50]
      ];
      for (var i = 0; i < asinh_of_exp_tests.length; i++) {
        (function() {
          var temp = asinh_of_exp_tests[i];
          u = temp[0];
          v = temp[1]
        })();
        y = math_util.asinh_of_exp(u);
        yy = v + u + Math.log(2.0);
        test.assert_rel_equal_eps(y, yy, (1.0e-16) * 5);
      };
      test.done(verbosity, "test_math_util");
