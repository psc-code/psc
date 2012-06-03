
#ifndef PARTICLES_CUDA_H
#define PARTICLES_CUDA_H

// FIXME, eventually this should probably be part of psc_particles_cuda_private.h

struct cuda_patch_flds {
  real *d_flds;
};

struct cuda_patch_prts {
  particles_cuda_dev_t d_part;
};

struct cuda_patch_ctx {
  struct cuda_patch_flds *h_cp_flds;
  struct cuda_patch_prts *h_cp_prts;
  int nr_patches;
  struct cuda_patch_flds *d_cp_flds;
  struct cuda_patch_prts *d_cp_prts;
  struct psc_particles **mprts_cuda;
  struct psc_fields **mflds_cuda;
};

EXTERN_C void cuda_patch_ctx_create(struct cuda_patch_ctx *cp, struct psc_mparticles *mprts,
				    struct psc_mfields *mflds);
EXTERN_C void cuda_patch_ctx_free(struct cuda_patch_ctx *cp);

#endif
