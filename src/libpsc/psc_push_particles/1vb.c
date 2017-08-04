
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
do_push_part_1vb_yz(struct psc_fields *flds, struct psc_mparticles *mprts, int p)
{
#ifdef PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  struct psc_particles *_prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(_prts);
#endif

  particle_range_t prts = particle_range_mprts(mprts, p);
  for (int b = 0; b < sub->nr_blocks; b++) {
    for (int n = sub->b_off[b]; n < sub->b_off[b+1]; n++) {
      push_one(prts.begin, n, flds, flds);
    }
  }
}

#else

static void
do_push_part_1vb_yz(struct psc_fields *flds, struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = particle_range_mprts(mprts, p);
  unsigned int n_prts = particle_range_size(prts);
  
  for (int n = 0; n < n_prts; n++) {
    push_one(prts.begin, n, flds, flds);
  }
}

#endif

void
SFX(psc_push_particles_push_mprts)(struct psc_push_particles *push,
				   struct psc_mparticles *mprts,
				   struct psc_mfields *mflds)
{
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_particles *_prts = psc_mparticles_get_patch(mprts, p);

    psc_fields_zero_range(flds, JXI, JXI + 3);
    ext_prepare_sort_before(_prts);
    do_push_part_1vb_yz(flds, mprts, p);
  }
}

