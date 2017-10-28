
#ifndef PARTICLES_CUDA_H
#define PARTICLES_CUDA_H

#define MAX_KINDS (4)

struct cuda_mfields_params {
  int mx[3];
  int ilg[3];
};

struct cuda_mparticles_params {
  float fnqs;
  float dxi[3];
  float b_dxi[3];
  int b_mx[3];
};

EXTERN_C void cuda_mfields_params_set(struct cuda_mfields_params *mflds_prm,
				      struct cuda_mfields *cmflds);
EXTERN_C void cuda_mparticles_params_set(struct cuda_mparticles_params *mprts_prm,
					 struct cuda_mparticles *cmprts, struct psc *psc);

// ======================================================================

EXTERN_C void cuda_mprts_find_block_indices_2(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_keys(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_indices_ids_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_scan_send_buf(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_bidx_to_key(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_check_ordered_total(struct psc_mparticles *mprts, int *n_prts_by_patch);

// FIXME, resolve this header mess eventually

#endif
