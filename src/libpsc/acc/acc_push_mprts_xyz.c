
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_XYZ
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define CALC_J CALC_J_1VB_SPLIT

#define SFX(s) s ## _xyz

#include "acc_push_mprts.c"

