
// ----------------------------------------------------------------------
// defaults

#ifndef CURR_CACHE
#define CURR_CACHE CURR_CACHE_NONE
#endif

typedef fields_t flds_curr_t;

// ----------------------------------------------------------------------

#define CURR_CACHE_GMEM 1
#define CURR_CACHE_N_REDUNDANT 1

#if CURR_CACHE == CURR_CACHE_NONE
#include "inc_curr_cache_none.c"
#elif CURR_CACHE == CURR_CACHE_SHIFT
#include "inc_curr_cache_shift.c"
#elif CURR_CACHE == CURR_CACHE_CUDA
#include "inc_curr_cache_cuda.c"
#elif CURR_CACHE == CURR_CACHE_CUDA2
#include "inc_curr_cache_cuda2.c"
#endif

// ----------------------------------------------------------------------

#if CALC_J == CALC_J_1VB_SPLIT
#include "inc_curr_1vb_split.c"
#elif CALC_J == CALC_J_1VB_VAR1
#include "inc_curr_1vb_var1.c"
#elif CALC_J == CALC_J_1VB_2D
#include "inc_curr_1vb_2d.c"
#endif

