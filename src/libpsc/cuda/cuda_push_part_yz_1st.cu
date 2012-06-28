
#include "psc_cuda.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define DIM DIM_YZ
#define PFX(x) yz4x4_1st_ ## x
#define CACHE_SHAPE_ARRAYS 7
#define CALC_CURRENT
#undef USE_SCRATCH
#define SW (2)

#include "cuda_common.h"
#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"
#include "common_curr.c"

#include "push_part_p1_1st.c"
#include "push_part_p2_1st.c"
#include "push_part_kernelcall_6.c"
