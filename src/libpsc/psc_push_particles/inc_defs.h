
// ======================================================================
// choices that determine which version of the pusher / deposition will
// be built

// ----------------------------------------------------------------------
// DIM

#define DIM_YZ 6
#define DIM_XYZ 7

// ----------------------------------------------------------------------
// CALC_J

#define CALC_J_1VB_SPLIT 1 // "original" V-B deposition with splitting along dims
#define CALC_J_1VB_VAR1  2 // V-B deposition variant with less code path divergence
#define CALC_J_1VB_2D    3 // V-B deposition variant with simpler out-of-plane current deposit

// ----------------------------------------------------------------------
// EM_CACHE

#define EM_CACHE_NONE 1
#define EM_CACHE_CUDA 2
#define EM_CACHE_CUDA2 3

// ----------------------------------------------------------------------
// CURR_CACHE

#define CURR_CACHE_NONE 1
#define CURR_CACHE_CUDA 2





