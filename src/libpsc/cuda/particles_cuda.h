
#ifndef PARTICLES_CUDA_H
#define PARTICLES_CUDA_H

// FIXME, eventually this should probably be part of psc_particles_cuda_private.h

struct cuda_patch {
  particles_cuda_dev_t d_part;
  real *d_flds;
};

__shared__ struct cuda_patch s_cpatch;

struct cuda_patch_ctx {
  struct cuda_patch *patch;
  int nr_patches;
  struct cuda_patch *d_patch;
  struct psc_particles **mprts_cuda;
  struct psc_fields **mflds_cuda;
};

EXTERN_C void cuda_patch_ctx_create(struct cuda_patch_ctx *cp, struct psc_mparticles *mprts,
				    struct psc_mfields *mflds);
EXTERN_C void cuda_patch_ctx_free(struct cuda_patch_ctx *cp);


#endif
