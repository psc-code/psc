
#include "psc.h"
#include "psc_particles_private.h"

#include <stdlib.h>
#include <assert.h>

const char *psc_topology_get_type(struct psc *psc, int p);

// ----------------------------------------------------------------------
// psc_mparticles_mix_setup

static void
psc_mparticles_mix_setup(struct psc_mparticles *mparticles)
{
  assert(mparticles->nr_particles_by_patch);

  mparticles->prts = calloc(mparticles->nr_patches, sizeof(*mparticles->prts));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    struct psc_particles *prts = psc_particles_create(MPI_COMM_SELF);
    psc_particles_set_type(prts, psc_topology_get_type(ppsc, p));
    char name[20]; sprintf(name, "prts%d", p);
    psc_particles_set_name(prts, name);
    prts->n_part = mparticles->nr_particles_by_patch[p];
    prts->flags = mparticles->flags;
    prts->p = p;
    psc_particles_setup(prts);
    mparticles->prts[p] = prts;
  }

  free(mparticles->nr_particles_by_patch);
  mparticles->nr_particles_by_patch = NULL;
}

// ======================================================================
// psc_mparticles: subclass "mix"
  
struct psc_mparticles_ops psc_mparticles_mix_ops = {
  .name                    = "mix",
  .setup                   = psc_mparticles_mix_setup,
};

