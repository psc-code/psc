
#ifndef PSC_PUSH_PARTICLES_PRIVATE_H
#define PSC_PUSH_PARTICLES_PRIVATE_H

#include <psc_push_particles.h>

struct psc_push_particles {
  struct mrc_obj obj;
};

struct psc_push_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_push_particles);
  void (*push_x)(struct psc_push_particles *push_particles,
		 mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_y)(struct psc_push_particles *push_particles,
		 mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_z)(struct psc_push_particles *push_particles,
		 mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_xy)(struct psc_push_particles *push_particles,
		  mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_xz)(struct psc_push_particles *push_particles,
		  mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_yz)(struct psc_push_particles *push_particles,
		  mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_xyz)(struct psc_push_particles *push_particles,
		   mparticles_base_t *particles, mfields_base_t *flds);

  void (*push_yz_a)(struct psc_push_particles *push_particles,
		    mparticles_base_t *particles, mfields_base_t *flds);
  void (*push_yz_b)(struct psc_push_particles *push_particles,
		    mparticles_base_t *particles, mfields_base_t *flds);
};

// ======================================================================

extern struct psc_push_particles_ops psc_push_particles_generic_c_ops;
extern struct psc_push_particles_ops psc_push_particles_1st_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_1sff_ops;
extern struct psc_push_particles_ops psc_push_particles_fortran_ops;
extern struct psc_push_particles_ops psc_push_particles_vay_ops;
extern struct psc_push_particles_ops psc_push_particles_sse2_ops;
extern struct psc_push_particles_ops psc_push_particles_cbe_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_1st_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_1vb_ops;

#define psc_push_particles_ops(push_particles) ((struct psc_push_particles_ops *)((push_particles)->obj.ops))

#endif
