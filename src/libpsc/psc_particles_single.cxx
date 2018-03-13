
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_particles_inc.h"
#include "psc_particles_double.h"

#if 0
static void _mrc_unused // FIXME
psc_particles_single_reorder(struct psc_mparticles *mprts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  if (!sub->need_reorder) {
    return;
  }
f
  int n_prts = patch->n_prts;
  for (int n = 0; n < n_prts; n++) {
    sub->particles_alt[n] = sub->particles[sub->b_ids[n]];
  }
  
  // swap in alt array
  particle_single_t *tmp = sub->particles;
  sub->particles = sub->particles_alt;
  sub->particles_alt = tmp;
  sub->need_reorder = false;
}
#endif

// ======================================================================
// psc_mparticles: subclass "single"

// ----------------------------------------------------------------------
// conversion to/from "double"

template<> const MparticlesBase::Map MparticlesSingle::convert_to_ = {
  { "copy_to_double"  , psc_mparticles_copy_to<MparticlesSingle, MparticlesDouble> },
};

template<> const MparticlesBase::Map MparticlesSingle::convert_from_ = {
  { "copy_from_double", psc_mparticles_copy_from<MparticlesSingle, MparticlesDouble> },
};

#include "psc_particles_common.cxx"

