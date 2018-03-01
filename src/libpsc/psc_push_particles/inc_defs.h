
#pragma once

#include "dim.hxx"

// ======================================================================
// choices that determine which version of the pusher / deposition will
// be built

// ----------------------------------------------------------------------
// DIM

#define DIM_X 1
#define DIM_Y 2
#define DIM_Z 4
#define DIM_XY (DIM_X | DIM_Y)
#define DIM_XZ (DIM_X | DIM_Z)
#define DIM_YZ (DIM_Y | DIM_Z)
#define DIM_XYZ (DIM_X | DIM_Y | DIM_Z)
#define DIM_1 0

// ----------------------------------------------------------------------
// CALC_J

#define CALC_J_1VB_SPLIT 1 // "original" V-B deposition with splitting along dims
#define CALC_J_1VB_VAR1  2 // V-B deposition variant with less code path divergence
#define CALC_J_1VB_2D    3 // V-B deposition variant with simpler out-of-plane current deposit

struct opt_calcj_esirkepov;
struct opt_calcj_1vb_split;
struct opt_calcj_1vb_var1;
struct opt_calcj_1vb_2d;

// ----------------------------------------------------------------------
// EM_CACHE

#define EM_CACHE_NONE 1
#define EM_CACHE_CUDA 2
#define EM_CACHE_CUDA2 3

// ----------------------------------------------------------------------
// ORDER

struct opt_order_1st {};
struct opt_order_2nd {};

struct opt_ip_1st;
struct opt_ip_1st_ec;
struct opt_ip_2nd;

// ----------------------------------------------------------------------


struct opt_ext_none;
struct opt_ext_prepare_sort;

