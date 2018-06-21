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
#endif

#if "LANG" eq "js"
#define IS_NAN(x) (isNaN(x))
#define IS_NONE(x) ((x)==null)
#define IS_REAL(x) (not (IS_NONE(x) or IS_NAN(x)))
#define NONE null
#define NAN (NaN)
# ... note that NAN==NAN is false
#define ARRAY_WITH_VAL(size,vvv) ((function() {a=new Array(size); for (var i=0; i<size; i++) {a[i]=vvv}; return a})())
# ... Surrounding this with parens makes it a function expression, which then gets evaluated on empty args ().
#     https://stackoverflow.com/a/9091416/1142217
#define EMPTY1D_BARE(size) ARRAY_WITH_VAL((size),0)
#define EMPTY2D_BARE(size) ARRAY_WITH_VAL((size),EMPTY1D_BARE((size)))
#define EMPTY3D_BARE(size) ARRAY_WITH_VAL((size),EMPTY2D_BARE((size)))
#define EMPTY1DIM(size) EMPTY1D_BARE(size)__NO_TRANSLATION__
#define EMPTY2DIM(size) EMPTY2D_BARE(size)__NO_TRANSLATION__
#define EMPTY3DIM(size) EMPTY3D_BARE(size)__NO_TRANSLATION__
#define PRINT print
# ... works in rhino and d8
#define IS_BROWSER (typeof window !== 'undefined')
# ... https://stackoverflow.com/q/26738943/1142217
#define THROW(message) throw message;
#endif

#endif
