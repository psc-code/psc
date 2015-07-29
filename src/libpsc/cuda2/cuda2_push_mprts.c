
#include "psc_cuda2.h"

// ----------------------------------------------------------------------
// FIXME
EXTERN_C void yz4x4_1vbec3d_gmem_cuda_push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds);

// end FIXME
// ----------------------------------------------------------------------

void
cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds,
			  struct psc_mparticles *mprts_cuda, struct psc_mfields *mflds_cuda)
{
  yz4x4_1vbec3d_gmem_cuda_push_mprts(mprts_cuda, mflds_cuda);
}

