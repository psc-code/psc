
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"
#include "1vb/psc_push_particles_1vb.h"

static void
psc_push_particles_1vbec_single_push_mprts(struct psc_push_particles *push,
					   struct psc_mparticles *mprts,
					   struct psc_mfields *mflds)
{
  int *gdims = ppsc->domain.gdims;

  if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
    psc_push_particles_1vbec_single_push_mprts_xyz(push, mprts, mflds);
  } else if (gdims[0] > 1 && gdims[1] == 1 && gdims[2] > 1) {
    psc_push_particles_1vbec_single_push_mprts_xyz_xz(push, mprts, mflds);
  } else if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
    psc_push_particles_1vbec_single_push_mprts_yz(push, mprts, mflds);
  } else if (gdims[0] == 1 && gdims[1] == 1 && gdims[2] == 1) {
    psc_push_particles_1vbec_single_push_mprts_1(push, mprts, mflds);
  } else {
    assert(0);
  }
}

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

struct psc_push_particles_ops psc_push_particles_1vbec_single_ops = {
  .name                  = "1vbec_single",
  .push_mprts            = psc_push_particles_1vbec_single_push_mprts,
  .stagger_mprts_yz      = psc_push_particles_1vbec_single_stagger_mprts_yz,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

