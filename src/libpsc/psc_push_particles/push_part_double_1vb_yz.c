
#include "psc_push_particles_double.h"

#define F3_CURR F3_C
#define F3_CACHE F3_C
#define F3_CACHE_TYPE "c"
#define psc_push_particles_push_yz psc_push_particles_1vb_double_push_a_yz

#define PUSHER_TYPE "1vb"
#define INTERPOLATE_1ST INTERPOLATE_1ST_STD

#include "1vb_yz.c"
