
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "inc_params.c"
#include "inc_cache.c"
#include "inc_interpolate.c"
#include "inc_push.c"
#include "inc_curr.c"
#include "inc_step.c"

// ======================================================================

#ifdef PUSHER_BY_BLOCK

static void
do_push_part_1vb_yz(struct psc_fields *flds, struct psc_particles *prts)
{
#ifdef PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);
#endif

  for (int b = 0; b < sub->nr_blocks; b++) {
    for (int n = sub->b_off[b]; n < sub->b_off[b+1]; n++) {
      push_one(prts, n, flds, flds);
    }
  }
}

#else

static void
do_push_part_1vb_yz(struct psc_fields *flds, struct psc_particles *prts)
{
  for (int n = 0; n < prts->n_part; n++) {
    push_one(prts, n, flds, flds);
  }
}

#endif

void
SFX(psc_push_particles_push_a)(struct psc_push_particles *push,
			       struct psc_particles *prts,
			       struct psc_fields *flds)
{
  psc_fields_zero_range(flds, JXI, JXI + 3);

  params_1vb_set(ppsc, NULL, NULL);
  ext_prepare_sort_before(prts);
  do_push_part_1vb_yz(flds, prts);
}

