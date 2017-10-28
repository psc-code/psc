
#ifndef PARTICLES_CUDA_H
#define PARTICLES_CUDA_H

#include "cuda_bits.h"

// ======================================================================

EXTERN_C void cuda_mprts_find_block_indices_2(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_keys(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_find_block_indices_ids_total(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_scan_send_buf(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_bidx_to_key(struct psc_mparticles *mprts);
EXTERN_C void cuda_mprts_check_ordered_total(struct psc_mparticles *mprts, int *n_prts_by_patch);

// FIXME, resolve this header mess eventually

#endif
