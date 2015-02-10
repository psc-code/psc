
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
  void (*push_mprts_yz)(struct psc_push_particles *push_particles,
			struct psc_mparticles *mprts, struct psc_mfields *mflds);
  
  unsigned int mp_flags; //< flags for _get_cuda(), alloc
};

// ======================================================================

extern struct psc_push_particles_ops psc_push_particles_generic_c_ops;
extern struct psc_push_particles_ops psc_push_particles_2nd_double_ops;
extern struct psc_push_particles_ops psc_push_particles_1p5_c_ops;
extern struct psc_push_particles_ops psc_push_particles_1st_ops;
extern struct psc_push_particles_ops psc_push_particles_1sff_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_mix_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_c_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_double_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ps_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ps2_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb2_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_double_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_double_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_single_by_block_ops;
extern struct psc_push_particles_ops psc_push_particles_fortran_ops;
extern struct psc_push_particles_ops psc_push_particles_vay_ops;
extern struct psc_push_particles_ops psc_push_particles_sse2_ops;
extern struct psc_push_particles_ops psc_push_particles_cbe_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_4x4_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_2x2_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_4x4_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_8x8_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_2x2_gmem_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_4x4_gmem_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_8x8_gmem_cuda_ops;

#define psc_push_particles_ops(push_particles) ((struct psc_push_particles_ops *)((push_particles)->obj.ops))

#endif
