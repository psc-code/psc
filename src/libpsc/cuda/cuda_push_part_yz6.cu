
#include "psc_cuda.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define DIM DIM_YZ
#define PFX(x) yz6_ ## x
#define CACHE_SHAPE_ARRAYS 5
#define CALC_CURRENT
#undef USE_SCRATCH
#define SW (2)

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"
#include "common_curr.c"

#include "push_part_p1.c"
#include "push_part_p2_noshift_6.c"
#include "push_part_kernelcall_6.c"
