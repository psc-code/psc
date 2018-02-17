
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

  for (int p = 0; p < mprts->n_patches(); p++) {
    mparticles_t::patch_t& prts = mprts[p];
    unsigned int n_prts = prts.size();
    for (int i = 0; i < n_prts; i++) {
      int j = random() % n_prts;
      std::swap(prts[i], prts[j]);
    }
  }

  mprts.put_as(mprts_base);
}

// ======================================================================
// psc_randomize: subclass "c"

struct psc_randomize_ops_double : psc_randomize_ops {
  psc_randomize_ops_double() {
    name                  = "c";
    run                   = psc_randomize_c_run;
  }
} psc_randomize_c_ops;
