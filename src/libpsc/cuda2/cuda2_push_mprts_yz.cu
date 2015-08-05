
#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_YZ
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#include "cuda2_push_mprts.cu"

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_yz

void
cuda2_1vbec_push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  int *bs = mprts_sub->bs;
  assert(bs[0] == BLOCKSIZE_X && bs[1] == BLOCKSIZE_Y && bs[2] == BLOCKSIZE_Z);
  cuda_push_mprts_ab(mprts, mflds);
}

