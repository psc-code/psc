
#ifndef PSC_CUDA2_H
#define PSC_CUDA2_H

#include "psc.h"

#ifdef __cplusplus

#define EXTERN_C extern "C"

#else

#define EXTERN_C

#endif

// ----------------------------------------------------------------------

static const int psc_particles_cuda2_bs[3] = { 4, 4, 4 };

EXTERN_C void *cuda_calloc(size_t nmemb, size_t size);
EXTERN_C void cuda_free(void *ptr);
EXTERN_C void cuda_memcpy_host_from_device(void *h_ptr, void *d_ptr, size_t n);
EXTERN_C void cuda_memcpy_device_from_host(void *d_ptr, void *h_ptr, size_t n);

EXTERN_C void psc_mfields_cuda2_copy_to_device(struct psc_mfields *mflds);
EXTERN_C void psc_mfields_cuda2_copy_to_host(struct psc_mfields *mflds);

EXTERN_C void psc_mparticles_cuda2_copy_to_device(struct psc_mparticles *mprts);
EXTERN_C void psc_mparticles_cuda2_copy_to_host(struct psc_mparticles *mprts);

EXTERN_C void cuda2_push_mflds_E_yz_gold(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_H_yz_gold(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_E_yz(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_H_yz(struct psc_mfields *mflds);

EXTERN_C void cuda2_push_mflds_E_xyz_gold(struct psc_mfields *mflds);
EXTERN_C void cuda2_push_mflds_H_xyz_gold(struct psc_mfields *mflds);

EXTERN_C void cuda2_1vbec_push_mprts_gold_yz(struct psc_mparticles *mprts,
					     struct psc_mfields *mflds);
EXTERN_C void cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts,
					struct psc_mfields *mflds);

EXTERN_C void cuda2_1vbec_push_mprts_gold_xyz(struct psc_mparticles *mprts,
					     struct psc_mfields *mflds);
EXTERN_C void cuda2_1vbec_push_mprts_xyz(struct psc_mparticles *mprts,
					 struct psc_mfields *mflds);

#endif


