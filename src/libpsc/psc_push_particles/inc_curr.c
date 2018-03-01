
#pragma once

// ----------------------------------------------------------------------

template<typename curr_cache_t, typename dim_t>
struct CurrentNone;

#include "inc_curr_1vb_split.c"
#include "inc_curr_1vb_var1.c"
#include "inc_curr_1vb_2d.c"

template<typename curr_cache_t, typename dim_t>
struct Current1vb;


