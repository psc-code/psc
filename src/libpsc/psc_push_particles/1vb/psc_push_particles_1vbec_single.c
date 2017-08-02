
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

// FIXME -> some header
void psc_push_particles_1vbec_single_push_mprts_yz(struct psc_push_particles *push,
						   struct psc_mparticles *mprts,
						   struct psc_mfields *mflds);
void psc_push_particles_1vbec_single_push_mprts_xyz(struct psc_push_particles *push,
						    struct psc_mparticles *mprts,
						    struct psc_mfields *mflds);

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

struct psc_push_particles_ops psc_push_particles_1vbec_single_ops = {
  .name                  = "1vbec_single",
  .push_mprts_yz         = psc_push_particles_1vbec_single_push_mprts_yz,
  .push_mprts_xyz        = psc_push_particles_1vbec_single_push_mprts_xyz,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

