
#include "psc.h"
#include "psc_particles_mix.h"

#include <stdlib.h>
#include <assert.h>

const char *psc_topology_get_type(struct psc *psc);

// ----------------------------------------------------------------------
// psc_mparticles_mix_create

static void
psc_mparticles_mix_create(struct psc_mparticles *mprts)
{
  struct psc_mparticles_mix *mix = psc_mparticles_mix(mprts);

  mix->sub = psc_mparticles_create(psc_mparticles_comm(mprts));
  psc_mparticles_add_child(mprts, (struct mrc_obj *) mix->sub);
  psc_mparticles_set_type(mix->sub, psc_topology_get_type(ppsc));
}

// ----------------------------------------------------------------------
// psc_mparticles_mix_setup

static void
psc_mparticles_mix_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_mix *mix = psc_mparticles_mix(mprts);

  // FIXME, this is just all way too hacky
  psc_mparticles_set_domain_nr_particles(mix->sub, mprts->domain,
					 mprts->nr_particles_by_patch);
  mprts->nr_particles_by_patch = NULL;
  psc_mparticles_setup_children(mprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_mix_setup_internals

static void
psc_mparticles_mix_setup_internals(struct psc_mparticles *mprts)
{
  struct psc_mparticles_mix *mix = psc_mparticles_mix(mprts);

  psc_mparticles_setup_internals(mix->sub);
}

// ======================================================================
// psc_mparticles: subclass "mix"
  
struct psc_mparticles_ops psc_mparticles_mix_ops = {
  .name                    = "mix",
  .size                    = sizeof(struct psc_mparticles_mix),
  .create                  = psc_mparticles_mix_create,
  .setup                   = psc_mparticles_mix_setup,
  .setup_internals         = psc_mparticles_mix_setup_internals,
};

