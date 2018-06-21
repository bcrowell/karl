#ifndef LANGUAGE_H
#define LANGUAGE_H

#if "LANG" eq "python"
#define IS_NAN(x) (math.isnan(x))
#define IS_NONE(x) ((x) is None)
#define IS_REAL(x) (not (IS_NONE(x) or IS_NAN(x)))
#define NONE (None)
#define NAN float('nan')
# ... note that NAN==NAN is false
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
#define EMPTY1DIM(x) [0,0,0,0,0]
# kludge: ignore x and allocate it with size 5, which is the biggest I ever need
#define EMPTY2DIM(x) [EMPTY1DIM(x),EMPTY1DIM(x),EMPTY1DIM(x),EMPTY1DIM(x),EMPTY1DIM(x)]
#define EMPTY3DIM(x) [EMPTY2DIM(x),EMPTY2DIM(x),EMPTY2DIM(x),EMPTY2DIM(x),EMPTY2DIM(x)]
#define PRINT print
# ... works in rhino and d8
#define IS_BROWSER (typeof window !== 'undefined')
# ... https://stackoverflow.com/q/26738943/1142217
#define THROW(message) throw message;
#endif

#endif
