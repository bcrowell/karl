#if "LANG" eq "python"
import test
from test import assert_rel_equal,assert_rel_equal_eps,assert_equal,assert_equal_eps,done,assert_boolean


verbosity=1 # default, can be modified

#endif

#if "LANG" eq "js"
load("test.js");
assert_boolean = test.assert_boolean;
assert_rel_equal = test.assert_rel_equal;
assert_equal = test.assert_equal;
assert_rel_equal_eps = test.assert_rel_equal_eps;
assert_equal_eps = test.assert_equal_eps;
#endif
