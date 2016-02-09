
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_YZ
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define CALC_J CALC_J_1VB_SPLIT

#define SFX(s) s ## _yz

#include "acc_push_mprts.c"

