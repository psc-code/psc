
#ifndef PSC_DEBUG_H
#define PSC_DEBUG_H

//#define DEBUG_FINITE

#ifdef DEBUG_FINITE
#define assert_finite(x) assert(isfinite(x))
#else
#define assert_finite(x) do {} while (0)
#endif

#endif
