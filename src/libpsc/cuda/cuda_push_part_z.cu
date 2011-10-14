
#include "psc_cuda.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define DIM DIM_Z
#define PFX(x) z_ ## x
#define CALC_CURRENT
#define SW (3)

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"
#include "common_curr.c"

#include "push_part_p1.c"
#include "push_part_p2.c"
#include "push_part_collect_currents.c"
#include "push_part_kernelcall.c"
