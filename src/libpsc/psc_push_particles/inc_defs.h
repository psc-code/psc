
#pragma once

#include "dim.hxx"

// ======================================================================
// choices that determine which version of the pusher / deposition will
// be built

// ----------------------------------------------------------------------
// CALC_J

#define CALC_J_1VB_SPLIT 1 // "original" V-B deposition with splitting along dims
#define CALC_J_1VB_VAR1  2 // V-B deposition variant with less code path divergence
#define CALC_J_1VB_2D    3 // V-B deposition variant with simpler out-of-plane current deposit

// ----------------------------------------------------------------------
// ORDER

struct opt_order_1st {};
struct opt_order_2nd {};

