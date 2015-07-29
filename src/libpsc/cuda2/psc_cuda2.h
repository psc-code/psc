
#ifndef PSC_CUDA2_H
#define PSC_CUDA2_H

#include "psc.h"

#ifdef __cplusplus

#define EXTERN_C extern "C"

#else

#define EXTERN_C

#endif


EXTERN_C void cuda2_push_mflds_E_yz(struct psc_mfields *mflds, struct psc_mfields *mflds_cuda);
EXTERN_C void cuda2_push_mflds_H_yz(struct psc_mfields *mflds);

EXTERN_C void cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts,
					struct psc_mfields *mflds);


static const int psc_particles_cuda2_bs[3] = { 1, 2, 2 };

#endif


