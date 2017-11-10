
#include "vpic_mparticles.h"

extern vpic_simulation *simulation;

// ======================================================================
// vpic_mparticles

// ----------------------------------------------------------------------
// vpic_mparticles_create

struct vpic_mparticles *
vpic_mparticles_create()
{
  return new vpic_mparticles;
}

// ----------------------------------------------------------------------
// vpic_mparticles_ctor_from_simulation

void
vpic_mparticles_ctor_from_simulation(struct vpic_mparticles *vmprts)
{
  vmprts->species_list = simulation->species_list;
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_nr_particles

int vpic_mparticles_get_nr_particles(struct vpic_mparticles *vmprts)
{
  int n_prts = 0;

  species_t *sp;
  LIST_FOR_EACH(sp, vmprts->species_list) {
    n_prts += sp->np;
  }

  return n_prts;
}

