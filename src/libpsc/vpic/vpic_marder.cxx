
#include "vpic_marder.h"

extern vpic_simulation *simulation;

// ======================================================================
// vpic_marder

// ----------------------------------------------------------------------
// vpic_marder_create

struct vpic_marder *
vpic_marder_create()
{
  return new vpic_marder;
}

// ----------------------------------------------------------------------
// vpic_marder_ctor_from_simulation

void
vpic_marder_ctor_from_simulation(struct vpic_marder *vmarder)
{
  vmarder->clean_div_e_interval = simulation->clean_div_e_interval;
  vmarder->clean_div_b_interval = simulation->clean_div_b_interval;
  vmarder->sync_shared_interval = simulation->sync_shared_interval;
  vmarder->num_div_e_round = simulation->num_div_e_round;
  vmarder->num_div_b_round = simulation->num_div_b_round;
}

// ----------------------------------------------------------------------
// vpic_marder_run

void
vpic_marder_run(struct vpic_marder *vmarder, struct vpic_mfields *vmflds,
		struct vpic_mparticles *vmprts, int step)
{
  // Divergence clean e
  int clean_div_e_interval = vmarder->clean_div_e_interval;
  if (clean_div_e_interval > 0 && step % clean_div_e_interval == 0) {
    vmarder->clean_div_e(vmflds, vmprts);
  }
  
  // Divergence clean b
  int clean_div_b_interval = vmarder->clean_div_b_interval;
  if (clean_div_b_interval > 0 && step % clean_div_b_interval == 0) {
    vmarder->clean_div_b(vmflds);
  }

  // Synchronize the shared faces
  int sync_shared_interval = vmarder->sync_shared_interval;
  if (sync_shared_interval > 0 && step % sync_shared_interval == 0) {
    vmarder->sync_faces(vmflds);
  }
}

