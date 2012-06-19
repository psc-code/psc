
#include "psc_push_particles_private.h"
#include "psc_push_particles_ps.h"

// FIXME, we should rather have our own calc_j
void psc_push_particles_1vb_single_calc_j_yz(struct psc_push_particles *push,
					     struct psc_particles *particles_base,
					     struct psc_fields *flds_base);

// ======================================================================
// psc_push_particles: subclass "1vb_ps"

struct psc_push_particles_ops psc_push_particles_1vb_ps_ops = {
  .name                  = "1vb_ps",
  .push_a_yz             = psc_push_particles_1vb_ps_push_a_yz,
  .calc_j_yz             = psc_push_particles_1vb_single_calc_j_yz,
};

// ======================================================================
// psc_push_particles: subclass "1vb_ps2"

struct psc_push_particles_ops psc_push_particles_1vb_ps2_ops = {
  .name                  = "1vb_ps2",
  .push_a_yz             = psc_push_particles_1vb_ps2_push_a_yz,
  .calc_j_yz             = psc_push_particles_1vb_single_calc_j_yz,
};

