
#pragma once

// ----------------------------------------------------------------------

template<typename _Order, typename _Dim, typename _fields_t>
struct CurrentNone;

#include "inc_curr_1vb_split.c"
#include "inc_curr_1vb_var1.c"
#include "inc_curr_1vb_2d.c"

template<typename _Order, typename _Dim, typename _fields_t>
struct Current1vb;


