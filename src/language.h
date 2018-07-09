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
#define APPEND_TO_ARRAY(x,y) ((x).append(y))
#define PRINT print
#define THROW(message) raise RuntimeError(message)
#define THROW_ARRAY(aaa) THROW((io_util.strcat(aaa)))
# ... usage: THROW_ARRAY(([...])) ... extra parens required by filepp so it believes it's a single argument
#define EXIT(err) os._exit(err)
#define BREAK break
#define MATH_E math.e
#define MATH_PI math.pi
#define CLONE_ARRAY_OF_FLOATS(x) copy.copy(x)
#          ... need extra parens in usages like CLONE_ARRAY_OF_FLOATS(([x,y])), otherwise filepp gets 
#              confused by comma and thinks this is multiple arguments
#define CLONE_ARRAY_OF_FLOATS2DIM(x) copy.copy(x)
#define HAS_KEY(x,y) ((y) in (x))
#define TRUE (True)
#define FALSE (False)
#define LEN(x) len(x)
#define CEIL(x) int(math.ceil(x))
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
#define APPEND_TO_ARRAY(x,y) ((x).push(y))
#define PRINT print
# ... works in rhino and d8
#define IS_BROWSER (typeof window !== 'undefined')
# ... https://stackoverflow.com/q/26738943/1142217
#define THROW(message) throw message;
#define THROW_ARRAY(aaa) THROW(io_util.strcat(aaa))
# ... usage: THROW_ARRAY(([...])) ... extra parens required by filepp so it believes it's a single argument
#define EXIT(err) (quit(err))
#define BREAK break
#define MATH_E Math.E
#define MATH_PI Math.PI
#define CLONE_ARRAY_OF_FLOATS(x) (karl.clone_array1d(x))
#           ... see notes above about usage with array literals
#define CLONE_ARRAY_OF_FLOATS2DIM(x) (karl.clone_array2d(x))
#define HAS_KEY(x,y) ((y) in (x))
#define TRUE (true)
#define FALSE (false)
#define LEN(x) ((x).length)
#define CEIL(x) math.ceil(x)
#                 ... works in rhino
#endif

#endif
