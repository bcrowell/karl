/* This library contains functions for cloning one-dimensional arrays,
   and for making new arrays initialized with zeroes. */

/* Clone a one-dimensional array. This is a shallow copy, so only use it
   on an array consisting of primitives such as floating-point numbers. */
karl.clone_array1d = function(a) {
  var n = a.length;
  var b = Array(n);
  for (var i=0; i<n; i++) {
    b[i] = a[i];
  }
  return b;
};

/* arrays initialized with zeroes */
karl.array1d = function(n) {
  a = Array(n);
  for (var i=0; i<n; i++) {
    a[i] = 0.0;
  }
  return a;
};
karl.array2d = function(m,n) {
  a = Array(m);
  for (var i=0; i<m; i++) {
    a[i] = [];
    for (j=0; j<n; j++) {
      a[i][j] = 0.0;
    }
  }
  return a;
};
karl.array3d = function(m,n,k) {
  a = Array(m);
  for (var i=0; i<m; i++) {
    a[i] = [];
    for (j=0; j<n; j++) {
      a[i][j] = [];
      for (k=0; k<n; k++) {
        a[i][j][k] = 0.0;
      }
    }
  }
  return a;
};
