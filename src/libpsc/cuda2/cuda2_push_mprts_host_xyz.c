
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_XYZ
#define BLOCKSIZE_X 4
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define EM_CACHE EM_CACHE_NONE
#define CALC_J CALC_J_1VB_SPLIT
#define F3_CURR F3_CUDA2

#define SFX(s) s ## _gold_xyz

#include "cuda2_push_mprts.cu"

