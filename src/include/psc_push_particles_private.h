
#ifndef PSC_PUSH_PARTICLES_PRIVATE_H
#define PSC_PUSH_PARTICLES_PRIVATE_H

#include <psc_push_particles.h>

struct psc_push_particles {
  struct mrc_obj obj;
};

struct psc_push_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_push_particles);
  void (*push_mprts)(struct psc_push_particles *push_particles,
		     struct psc_mparticles *mprts, struct psc_mfields *mflds);
  void (*stagger_mprts)(struct psc_push_particles *push_particles,
			struct psc_mparticles *mprts, struct psc_mfields *mflds);
  void (*prep)(struct psc_push_particles *push_particles,
	       struct psc_mparticles *mprts, struct psc_mfields *mflds);
  void (*stagger_mprts_yz)(struct psc_push_particles *push_particles,
			   struct psc_mparticles *mprts, struct psc_mfields *mflds);
  void (*stagger_mprts_1)(struct psc_push_particles *push_particles,
			  struct psc_mparticles *mprts, struct psc_mfields *mflds);
  
  unsigned int mp_flags; //< flags for _get_as(CUDA), alloc
  const char *particles_type;
};

// ======================================================================

#define psc_push_particles_ops(push_particles) ((struct psc_push_particles_ops *)((push_particles)->obj.ops))

#endif
