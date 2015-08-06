
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


// ======================================================================

#ifdef __CUDACC__
#define CUDA_DEVICE __device__
#define CUDA_CONSTANT __constant__ __device__
#else
#define CUDA_DEVICE
#define CUDA_CONSTANT
#endif
