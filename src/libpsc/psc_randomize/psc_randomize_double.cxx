
#include "psc_randomize_private.h"
#include "psc_particles_as_double.h"

#include <mrc_profile.h>

// FIXME, this is very inefficient if followed by a sort..

// ----------------------------------------------------------------------
// psc_randomize_c_run

static void
psc_randomize_c_run(struct psc_randomize *randomize,
		    struct psc_mparticles *mprts_base)
{
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();

  for (int p = 0; p < mprts.n_patches(); p++) {
    particle_range_t prts = particle_range_mprts(mprts.mprts(), p);
    unsigned int n_prts = particle_range_size(prts);
    for (int i = 0; i < n_prts; i++) {
      int j = random() % n_prts;
      particle_t tmp = *particle_iter_at(prts.begin, i);
      *particle_iter_at(prts.begin, i) = *particle_iter_at(prts.begin, j);
      *particle_iter_at(prts.begin, j) = tmp;
    }
  }

  mprts.put_as(mprts_base);
}

// ======================================================================
// psc_randomize: subclass "c"

struct psc_randomize_ops psc_randomize_c_ops = {
  .name                  = "c",
  .run                   = psc_randomize_c_run,
};
