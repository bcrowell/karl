#ifndef LANGUAGE_H
#define LANGUAGE_H

#if "LANG" eq "python"
#define __NO_TRANSLATION__
#...a no-op, only has significance if we're translating to js
#define IS_NAN(x) (math.isnan(x))
#define IS_NONE(x) ((x) is None)
#define IS_REAL(x) (not (IS_NONE(x) or IS_NAN(x)))
#define NONE (None)
#define NAN float('nan')
# ... note that NAN==NAN is false
#define EMPTY1DIM(x) [0 for i in range(x)]
#define EMPTY2DIM(x) [[0 for i in range(x)] for j in range(x)]
#define EMPTY3DIM(x) [[[0 for i in range(x)] for j in range(x)] for k in range(x)]
#define PRINT print
#define THROW(message) raise RuntimeError(message)
#define MATH_E math.e
#endif

#if "LANG" eq "js"
#define IS_NAN(x) (isNaN(x))
#define NAN (NaN)
# ... note that NAN==NAN is false, so use IS_NAN
#define IS_NONE(x) ((x)==null)
#define IS_REAL(x) (not (IS_NONE(x) or IS_NAN(x)))
#define NONE null
karl.load("lib/array");
#define EMPTY1DIM(x) karl.array1d((x))
#define EMPTY2DIM(x) karl.array2d((x),(x))
#define EMPTY3DIM(x) karl.array3d((x),(x),(x))
#define PRINT print
# ... works in rhino and d8
#define IS_BROWSER (typeof window !== 'undefined')
# ... https://stackoverflow.com/q/26738943/1142217
#define THROW(message) throw message;
#define MATH_E Math.E
#endif

#endif
