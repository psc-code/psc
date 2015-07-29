
#ifndef PSC_CUDA2_H
#define PSC_CUDA2_H

#include "psc.h"

#ifdef __cplusplus

#define EXTERN_C extern "C"

#else

#define EXTERN_C

#endif

static const int psc_particles_cuda2_bs[3] = { 1, 2, 2 };

EXTERN_C void *cuda_calloc(size_t nmemb, size_t size);
EXTERN_C void cuda_free(void *ptr);
EXTERN_C void cuda_memcpy_host_from_device(void *h_ptr, void *d_ptr, size_t n);
EXTERN_C void cuda_memcpy_device_from_host(void *d_ptr, void *h_ptr, size_t n);

EXTERN_C void cuda2_push_mflds_E_yz_gold(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_H_yz_gold(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_E_yz(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_H_yz(struct psc_mfields *mflds);

EXTERN_C void cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts,
					struct psc_mfields *mflds);

#endif


