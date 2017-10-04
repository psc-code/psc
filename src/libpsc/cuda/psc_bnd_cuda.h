
#ifndef PSC_BND_CUDA_H
#define PSC_BND_CUDA_H

#include "psc_bnd_private.h"
#include "psc_particles_cuda.h"

#include "cuda_wrap.h"

/* EXTERN_C void cuda_copy_offsets_to_dev(struct psc_particles *prts, unsigned int *h_offsets); */
/* EXTERN_C void cuda_find_block_indices(struct psc_particles *prts, unsigned int *d_bidx); */
/* EXTERN_C void cuda_find_block_indices_2(struct psc_particles *prts, unsigned int *d_bidx, */
/* 					int start); */
/* EXTERN_C void cuda_find_block_indices_3(struct psc_particles *prts, unsigned int *d_bidx, */
/* 					unsigned int *d_alt_bidx, */
/* 					int start, unsigned int *bn_idx, */
/* 					unsigned int *bn_off); */
/* EXTERN_C int cuda_exclusive_scan_2(struct psc_particles *prts, unsigned int *_d_vals, */
/* 				   unsigned int *_d_sums); */
/* EXTERN_C void cuda_reorder_send_buf(int p, struct psc_particles *prts,  */
/* 				    unsigned int *d_bidx, unsigned int *d_sums, */
/* 				    int n_send); */
/* EXTERN_C void cuda_reorder(struct psc_particles *prts, unsigned int *d_ids); */


// FIXME separate header

EXTERN_C void *sort_pairs_create(const int b_mx[3]);
EXTERN_C void sort_pairs_destroy(void *_sp);
EXTERN_C void sort_pairs_device_2(void *_sp, unsigned int *d_bidx,
				  unsigned int *d_alt_ids,
				  int n_part, int *d_offsets);

EXTERN_C void *sort_pairs_3_create(const int b_mx[3]);
EXTERN_C void sort_pairs_3_destroy(void *_sp);
EXTERN_C void sort_pairs_3_device(void *_sp, unsigned int *d_bidx,
				  unsigned int *d_alt_bidx, unsigned int *d_alt_ids,
				  int n_part, int *d_offsets,
				  int n_part_prev, unsigned int *bn_cnts);

#endif

