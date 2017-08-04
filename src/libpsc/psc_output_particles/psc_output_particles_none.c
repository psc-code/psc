
#include "psc_output_particles_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_output_particles_none_run

static void
psc_output_particles_none_run(struct psc_output_particles *out,
				 struct psc_mparticles *particles_base)
{
}

// ======================================================================
// psc_output_particles: subclass "none"

struct psc_output_particles_ops psc_output_particles_none_ops = {
  .name                  = "none",
  .run                   = psc_output_particles_none_run,
};
