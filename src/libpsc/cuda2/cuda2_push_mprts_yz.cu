
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_YZ

#include "cuda2_push_mprts.cu"

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_yz

void
cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  int *bs = mprts_sub->bs;
  if (bs[0] == 1 && bs[1] == 4 && bs[2] == 4) {
    cuda_push_mprts_ab<1, 4, 4>(mprts, mflds);
  } else {
    mprintf("unknown bs %d %d %d\n", bs[0], bs[1], bs[2]);
  }
}

