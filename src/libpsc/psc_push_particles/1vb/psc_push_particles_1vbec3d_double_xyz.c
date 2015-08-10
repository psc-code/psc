
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "../inc_defs.h"

#define DIM DIM_XYZ
#define CALC_J CALC_J_1VB_SPLIT
#define F3_CURR F3_C
//#define F3_CACHE F3_C
//#define F3_CACHE_TYPE "c"
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC

#define psc_push_particles_push_a_xyz psc_push_particles_1vbec3d_double_push_a_xyz

#include "1vb_xyz.c"

