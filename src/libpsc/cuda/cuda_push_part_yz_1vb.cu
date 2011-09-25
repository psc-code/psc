
#include "psc_cuda.h"

#define DIM DIM_YZ
#define PFX(x) yz_1vb_ ## x
#define CACHE_SHAPE_ARRAYS 7
#define CALC_CURRENT
#undef USE_SCRATCH
#define SW (2)

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"
#include "common_curr.c"

#include "push_part_p1_1st.c"
#include "push_part_p2_1vb.c"
#include "push_part_kernelcall_6.c"
