
#include "psc_cuda.h"

#define DIM DIM_YZ
#define PFX(x) yz4_ ## x
#define CALC_CURRENT
#define SW (2)

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"
#include "common_curr.c"

#include "push_part_p1.c"
#include "push_part_p2_noshift_2.c"
#include "push_part_kernelcall.c"
