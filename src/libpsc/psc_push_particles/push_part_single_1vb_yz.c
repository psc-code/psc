
#include "psc_push_particles_single.h"

#include "psc_fields_single.h"

#define F3_CURR F3_S
#define F3_CACHE F3_S
#define psc_push_particles_1vb_push_yz psc_push_particles_single_1vb_push_yz

#include "1vb_yz.c"
