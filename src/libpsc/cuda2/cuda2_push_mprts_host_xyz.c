
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_XYZ
#define BLOCKSIZE_X 8
#define BLOCKSIZE_Y 8
#define BLOCKSIZE_Z 8
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define EM_CACHE EM_CACHE_CUDA
#define CURR_CACHE CURR_CACHE_CUDA
#define CALC_J CALC_J_1VB_SPLIT

#define SFX(s) s ## _gold_xyz

#include "cuda2_push_mprts.cu"

