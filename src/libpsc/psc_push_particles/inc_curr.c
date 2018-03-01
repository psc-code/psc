
// ----------------------------------------------------------------------

template<typename curr_cache_t, typename dim_t>
struct CurrentNone;

#if CALC_J == CALC_J_1VB_SPLIT
#include "inc_curr_1vb_split.c"
#elif CALC_J == CALC_J_1VB_VAR1
#include "inc_curr_1vb_var1.c"
#elif CALC_J == CALC_J_1VB_2D
#include "inc_curr_1vb_2d.c"
#else
template<typename curr_cache_t, typename dim_t>
struct Current1vb;
#endif

