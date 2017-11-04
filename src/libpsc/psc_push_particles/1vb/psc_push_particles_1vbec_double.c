
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

// FIXME -> some header
void psc_push_particles_1vbec_double_push_mprts_yz(struct psc_push_particles *push,
						   struct psc_mparticles *mprts,
						   struct psc_mfields *mflds);
void psc_push_particles_1vbec_double_push_mprts_xyz(struct psc_push_particles *push,
						    struct psc_mparticles *mprts,
						    struct psc_mfields *mflds);
void psc_push_particles_1vbec_double_push_mprts_1(struct psc_push_particles *push,
						  struct psc_mparticles *mprts,
						  struct psc_mfields *mflds);
void psc_push_particles_1vbec_double_stagger_mprts_1(struct psc_push_particles *push,
						     struct psc_mparticles *mprts,
						     struct psc_mfields *mflds);

static void
psc_push_particles_1vbec_double_push_mprts(struct psc_push_particles *push,
					   struct psc_mparticles *mprts,
					   struct psc_mfields *mflds)
{
  int *gdims = ppsc->domain.gdims;

  if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
    psc_push_particles_1vbec_double_push_mprts_xyz(push, mprts, mflds);
  } else if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
    psc_push_particles_1vbec_double_push_mprts_yz(push, mprts, mflds);
  } else if (gdims[0] == 1 && gdims[1] == 1 && gdims[2] == 1) {
    psc_push_particles_1vbec_double_push_mprts_1(push, mprts, mflds);
  } else {
    assert(0);
  }
}

// ======================================================================
// psc_push_particles: subclass "1vbec_double"

struct psc_push_particles_ops psc_push_particles_1vbec_double_ops = {
  .name                  = "1vbec_double",
  .push_mprts            = psc_push_particles_1vbec_double_push_mprts,
  //  .stagger_mprts_1       = psc_push_particles_1vbec_double_stagger_mprts_1,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

