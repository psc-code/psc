
#include "psc_randomize_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_randomize_fortran_run

static void
psc_randomize_fortran_run(struct psc_randomize *randomize,
			  mparticles_base_t *particles_base)
{
  assert(psc.nr_patches == 1);
  
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_randomize", 1., 0, 0);
  }
  prof_start(pr);
  PIC_randomize(&particles.p[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
}

// ======================================================================
// psc_randomize: subclass "fortran"

struct psc_randomize_ops psc_randomize_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_randomize_fortran_run,
};
