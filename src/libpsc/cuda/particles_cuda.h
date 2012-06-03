
#ifndef PARTICLES_CUDA_H
#define PARTICLES_CUDA_H

// FIXME, eventually this should probably be part of psc_particles_cuda_private.h

struct cuda_patch_flds {
  real *d_flds;
};

struct cuda_patch_prts {
  particles_cuda_dev_t d_part;
};

struct cuda_mprts {
  struct cuda_patch_prts *h_cp_prts;
  struct cuda_patch_prts *d_cp_prts;
  int nr_patches;
  struct psc_particles **mprts_cuda;
};

struct cuda_mflds {
  struct cuda_patch_flds *h_cp_flds;
  struct cuda_patch_flds *d_cp_flds;
  int nr_patches;
  struct psc_fields **mflds_cuda;
};

EXTERN_C void cuda_mprts_create(struct cuda_mprts *cuda_mprts, struct psc_mparticles *mprts);
EXTERN_C void cuda_mflds_create(struct cuda_mflds *cuda_mflds, struct psc_mfields *mflds);
EXTERN_C void cuda_mprts_destroy(struct cuda_mprts *cuda_mprts);
EXTERN_C void cuda_mflds_destroy(struct cuda_mflds *cuda_mflds);

#define MAX_KINDS (4)

struct cuda_params {
  real dt;
  real dxi[3];
  real dqs;
  real fnqs;
  real fnqys, fnqzs;
  int mx[3];
  int ilg[3];
  int b_mx[3];
  int *d_error_count;
  real dq[MAX_KINDS];
};

EXTERN_C void set_params(struct cuda_params *prm, struct psc *psc,
			 struct psc_particles *prts, struct psc_fields *pf);
EXTERN_C void free_params(struct cuda_params *prm);

#endif
