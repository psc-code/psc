
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define F3_CURR F3_S
#define F3_CACHE F3_S
#define F3_CACHE_TYPE "single"

#define PUSHER_TYPE "1vbec3d"
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC

#include "1vb_yz.c"

// ======================================================================
// psc_push_particles: subclass "1vbec3d_single"

struct psc_push_particles_ops psc_push_particles_1vbec3d_single_ops = {
  .name                  = "1vbec3d_single",
  .push_a_yz             = psc_push_particles_push_a_yz,
};

