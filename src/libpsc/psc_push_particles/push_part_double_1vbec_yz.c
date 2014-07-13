
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#define F3_CURR F3_C
#define F3_CACHE F3_C
#define F3_CACHE_TYPE "c"

#define PUSHER_TYPE "1vbec"
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define VB_2D

#include "1vb_yz.c"

// ======================================================================
// psc_push_particles: subclass "1vbec_double"

struct psc_push_particles_ops psc_push_particles_1vbec_double_ops = {
  .name                  = "1vbec_double",
  .push_a_yz             = psc_push_particles_push_a_yz,
};

