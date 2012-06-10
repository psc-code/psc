
#ifndef PARTICLES_CUDA_H
#define PARTICLES_CUDA_H

// FIXME, eventually this should probably be part of psc_particles_cuda_private.h

struct cuda_patch_flds {
  real *d_flds;
};

struct cuda_mflds {
  struct cuda_patch_flds *h_cp_flds;
  struct cuda_patch_flds *d_cp_flds;
  int nr_patches;
  struct psc_fields **mflds_cuda;
};

EXTERN_C void psc_mparticles_cuda_copy_to_dev(struct psc_mparticles *mprts);
EXTERN_C void cuda_mflds_create(struct cuda_mflds *cuda_mflds, struct psc_mfields *mflds);
EXTERN_C void cuda_mflds_destroy(struct cuda_mflds *cuda_mflds);

#define MAX_KINDS (4)

struct cuda_params {
  real dt;
  real dxi[3];
  real b_dxi[3];
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

// ======================================================================

EXTERN_C void cuda_mprts_sort_initial(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_reorder_and_offsets(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_indices_2(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_indices_2_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_indices_ids_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_scan_send_buf(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_spine_reduce(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_n_send(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_scan_send_buf_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_reorder_send_buf_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_copy_from_dev(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_convert_from_cuda(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_copy_to_dev(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_indices_3(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_sort(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_reorder(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_check_ordered_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_sort_pairs_device(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_do_sort(struct psc_mparticles *mprts);

// FIXME, resolve this header mess eventually

EXTERN_C void __psc_mparticles_cuda_setup(struct psc_mparticles *mprts);
EXTERN_C void __psc_mparticles_cuda_free(struct psc_mparticles *mprts);

#endif
