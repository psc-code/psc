
#include "psc_cuda2.h"

#define DIM DIM_XYZ
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define CALC_J CALC_J_1VB_SPLIT
#define F3_CURR F3_CUDA2

#define SFX(s) s ## _xyz

#include "cuda2_push_mprts_host.c"

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts_gold

void
SFX(cuda2_1vbec_push_mprts_gold)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  params_1vb_set(ppsc, 0);

  psc_mfields_zero_range(mflds, JXI, JXI + 3);

  for (int b = 0; b < mprts_sub->nr_blocks_total; b++) {
    int p = b / mprts_sub->nr_blocks;
    for (int n = mprts_sub->h_b_off[b]; n < mprts_sub->h_b_off[b+1]; n++) {
      push_one(mprts, mflds, n, p);
    }
  }
}

