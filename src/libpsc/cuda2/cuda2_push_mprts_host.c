
#include "psc_cuda2.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define F3_CURR F3_S
#define F3_CACHE F3_S
#define F3_CACHE_TYPE "single"

#define INTERPOLATE_1ST INTERPOLATE_1ST_EC

#include "../psc_push_particles/1vb_yz.c"

void
cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    psc_push_particles_push_a_yz(NULL, psc_mparticles_get_patch(mprts, p),
				 psc_mfields_get_patch(mflds, p));
  }
}

