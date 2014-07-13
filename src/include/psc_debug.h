
#ifndef PSC_DEBUG_H
#define PSC_DEBUG_H

//#define DEBUG_FINITE

#ifdef DEBUG_FINITE
#define assert_finite(x) assert(isfinite(x))
#else
#define assert_finite(x) do {} while (0)
#endif

// ======================================================================
// for NaN poisoning

static inline void
float_set_nan(float *x)
{
  *(int *) x = 0x7fbfffff;
}

static inline void
double_set_nan(double *x)
{
  *(long long *) x = 0x7ff8000000000000;
}

#endif
