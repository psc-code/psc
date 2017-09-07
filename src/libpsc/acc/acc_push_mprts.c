
#include "psc_acc.h"
#include "psc_particles_as_acc.h"
#include "psc_fields_as_acc.h"

typedef struct psc_fields *curr_cache_t;
#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_cache.c"
#include "../psc_push_particles/inc_interpolate.c"
#include "../psc_push_particles/inc_push.c"
#include "../psc_push_particles/inc_curr.c"
#include "../psc_push_particles/inc_step.c"

// ----------------------------------------------------------------------
// push_mprts_loop

static void
push_mprts_loop(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_acc *mprts_sub = psc_mparticles_acc(mprts);
  mprts_array_t mprts_arr = { .xi4 = mprts_sub->xi4, .pxi4 = mprts_sub->pxi4, };

  for (int b = 0; b < mprts_sub->nr_blocks_total; b++) {
    int p = b / mprts_sub->nr_blocks;
    for (int n = mprts_sub->b_off[b]; n < mprts_sub->b_off[b+1]; n++) {
     struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
     push_one(mprts_arr, n, flds, flds);
    }
  }
}

// ----------------------------------------------------------------------
// acc_1vbec_push_mprts

void
SFX(acc_1vbec_push_mprts)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  c_prm_set(ppsc);
  params_1vb_set(ppsc, mprts, mflds);
  psc_mfields_zero_range(mflds, JXI, JXI + 3);
  push_mprts_loop(mprts, mflds);
}

