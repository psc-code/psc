
#include "psc_push_particles_double.h"

#define F3_CURR F3_C
#define F3_CACHE F3_C
#define cache_fields_from_em cache_fields_c_from_em
#define cache_fields_to_j cache_fields_c_to_j
#define psc_push_particles_1vb_push_yz psc_push_particles_double_1vb_push_yz

#include "1vb_yz.c"
