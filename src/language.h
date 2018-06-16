#if "LANG" eq "python"
#define IS_NAN(x) math.isnan(x)
#define IS_NONE(x) ((x) is not None)
#define IS_REAL(x) (not (IS_NONE(x) or IS_NAN(x)))
#define NAN float('nan')
  # ... note that NAN==NAN is false
#endif
