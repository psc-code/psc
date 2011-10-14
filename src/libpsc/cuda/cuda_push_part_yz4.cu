
#include "psc_cuda.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define DIM DIM_YZ
#define PFX(x) yz4_ ## x
#define CACHE_SHAPE_ARRAYS 4
#define CALC_CURRENT
#define USE_SCRATCH
#define SW (2)

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"
#include "common_curr.c"

#include "push_part_p1.c"
#include "push_part_p2_noshift_2.c"
#include "push_part_collect_currents.c"
#include "push_part_kernelcall.c"
