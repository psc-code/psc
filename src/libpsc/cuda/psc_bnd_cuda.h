
#ifndef PSC_BND_CUDA_H
#define PSC_BND_CUDA_H

#include "psc_bnd_private.h"
#include "psc_particles_cuda.h"

#include "cuda_bits.h"

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

