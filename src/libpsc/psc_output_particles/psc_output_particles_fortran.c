
#include "psc_output_particles_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_output_particles_fortran_run

static void
psc_output_particles_fortran_run(struct psc_output_particles *out,
				 mparticles_base_t *particles_base)
{
  assert(ppsc->nr_patches == 1);
  static int pr;
  if (!pr) {
    pr = prof_register("fort_out_part", 1., 0, 0);
  }
  prof_start(pr);
  mparticles_fortran_t particles;
  mparticles_fortran_get(&particles, particles_base);

  OUT_part(&particles.p[0]);
  
  mparticles_fortran_put(&particles, particles_base);
  prof_stop(pr);
}

// ======================================================================
// psc_output_particles: subclass "fortran"

struct psc_output_particles_ops psc_output_particles_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_output_particles_fortran_run,
};
