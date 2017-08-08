
#include "psc.h"
#include "psc_particles_private.h"

#include <mrc_profile.h>
#include <string.h>

// ======================================================================
// psc_particles class

struct mrc_class_psc_particles mrc_class_psc_particles = {
  .name             = "psc_particles",
  .size             = sizeof(struct psc_particles),
};

