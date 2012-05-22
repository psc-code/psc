
#ifndef PSC_PUSH_PARTICLES_PRIVATE_H
#define PSC_PUSH_PARTICLES_PRIVATE_H

#include <psc_push_particles.h>

struct psc_push_particles {
  struct mrc_obj obj;
};

struct psc_push_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_push_particles);
  void (*push_a_x)(struct psc_push_particles *push_particles,
		   struct psc_particles *prts, struct psc_fields *flds);
  void (*push_a_y)(struct psc_push_particles *push_particles,
		   struct psc_particles *prts, struct psc_fields *flds);
  void (*push_a_z)(struct psc_push_particles *push_particles,
		   struct psc_particles *prts, struct psc_fields *flds);
  void (*push_a_xy)(struct psc_push_particles *push_particles,
		    struct psc_particles *prts, struct psc_fields *flds);
  void (*push_a_xz)(struct psc_push_particles *push_particles,
		    struct psc_particles *prts, struct psc_fields *flds);
  void (*push_a_yz)(struct psc_push_particles *push_particles,
		    struct psc_particles *prts, struct psc_fields *flds);
  void (*push_b_yz)(struct psc_push_particles *push_particles,
		    struct psc_particles *prts, struct psc_fields *flds);
  void (*push_a_xyz)(struct psc_push_particles *push_particles,
		     struct psc_particles *prts, struct psc_fields *flds);

  void (*calc_j_x)(struct psc_push_particles *push_particles,
		   mparticles_base_t *particles, mfields_base_t *flds);
  void (*calc_j_y)(struct psc_push_particles *push_particles,
		   mparticles_base_t *particles, mfields_base_t *flds);
  void (*calc_j_z)(struct psc_push_particles *push_particles,
		   mparticles_base_t *particles, mfields_base_t *flds);
  void (*calc_j_xy)(struct psc_push_particles *push_particles,
		    mparticles_base_t *particles, mfields_base_t *flds);
  void (*calc_j_xz)(struct psc_push_particles *push_particles,
		    mparticles_base_t *particles, mfields_base_t *flds);
  void (*calc_j_yz)(struct psc_push_particles *push_particles,
		    mparticles_base_t *particles, mfields_base_t *flds);
  void (*calc_j_xyz)(struct psc_push_particles *push_particles,
		     mparticles_base_t *particles, mfields_base_t *flds);
  
  unsigned int mp_flags; //< flags for _get_cuda(), alloc
};

// ======================================================================

extern struct psc_push_particles_ops psc_push_particles_generic_c_ops;
extern struct psc_push_particles_ops psc_push_particles_1st_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_1sff_ops;
extern struct psc_push_particles_ops psc_push_particles_single_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_double_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_ps_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_fortran_ops;
extern struct psc_push_particles_ops psc_push_particles_vay_ops;
extern struct psc_push_particles_ops psc_push_particles_sse2_ops;
extern struct psc_push_particles_ops psc_push_particles_cbe_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_1st_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_2x2_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_cuda_8x8_1vb_ops;

#define psc_push_particles_ops(push_particles) ((struct psc_push_particles_ops *)((push_particles)->obj.ops))

#endif
