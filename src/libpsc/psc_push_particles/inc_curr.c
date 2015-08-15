
#define CURR_CACHE_GMEM 1
#define CURR_CACHE_N_REDUNDANT 1

#if PSC_FIELDS_AS_CUDA2

#if CURR_CACHE == CURR_CACHE_NONE
#include "inc_curr_cache_none.c"
#elif CURR_CACHE == CURR_CACHE_CUDA
#include "inc_curr_cache_cuda.c"
#elif CURR_CACHE == CURR_CACHE_CUDA2
#include "inc_curr_cache_cuda2.c"
#endif

#else

typedef struct psc_fields * curr_cache_t;

CUDA_DEVICE static inline void
curr_add(curr_cache_t curr_cache, int m, int jx, int jy, int jz, real val)
{
  F3_CURR(curr_cache, JXI+m, jx,jy,jz) += val;
}

#endif

// ----------------------------------------------------------------------

#if CALC_J == CALC_J_1VB_SPLIT
#include "inc_curr_1vb_split.c"
#elif CALC_J == CALC_J_1VB_VAR1
#include "inc_curr_1vb_var1.c"
#elif CALC_J == CALC_J_1VB_2D
#include "inc_curr_1vb_2d.c"
#endif

