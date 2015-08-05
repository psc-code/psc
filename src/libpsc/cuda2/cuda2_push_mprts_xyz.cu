
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_XYZ

#include "cuda2_push_mprts.cu"

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_xyz

void
cuda2_1vbec_push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  int *bs = mprts_sub->bs;
  if (bs[0] == 4 && bs[1] == 4 && bs[2] == 4) {
    cuda_push_mprts_ab<4, 4, 4>(mprts, mflds);
  } else {
    mprintf("unknown bs %d %d %d\n", bs[0], bs[1], bs[2]);
  }
}


