
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "inc_params.c"
#include "inc_interpolate.c"
#include "inc_push.c"
#include "inc_curr.c"
#include "c_common_push.c"

// ======================================================================
// EXT_PREPARE_SORT
//
// if enabled, calculate the new block position for each particle once
// it's moved, and also append particles that left the block to an extra
// list at the end of all local particles (hopefully there's enough room...)

#ifdef EXT_PREPARE_SORT

#ifdef PSC_PARTICLES_AS_SINGLE
static struct psc_particles_single *prts_sub;
#pragma omp threadprivate(prts_sub)
#endif

static inline void
ext_prepare_sort_before(struct psc_particles *prts)
{
  prts_sub = psc_particles_single(prts);
  memset(prts_sub->b_cnt, 0,
	 (prts_sub->nr_blocks + 1) * sizeof(*prts_sub->b_cnt));
}

static inline void
ext_prepare_sort(struct psc_particles *prts, int n, particle_t *prt,
		 int *b_pos)
{
  /* FIXME, only if blocksize == 1! */
  int *b_mx = prts_sub->b_mx;
  if (b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
      b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
    prts_sub->b_idx[n] = b_pos[2] * b_mx[1] + b_pos[1];
  } else { /* out of bounds */
    prts_sub->b_idx[n] = prts_sub->nr_blocks;
    assert(prts_sub->b_cnt[prts_sub->nr_blocks] < prts_sub->n_alloced);
    /* append to back */
    *particles_get_one(prts, prts->n_part + prts_sub->b_cnt[prts_sub->nr_blocks]) = *prt;
  }
  prts_sub->b_cnt[prts_sub->b_idx[n]]++;
}

#else

static inline void
ext_prepare_sort_before(struct psc_particles *prts)
{
}

static inline void
ext_prepare_sort(struct psc_particles *prts, int n, particle_t *prt,
		 int *b_pos)
{
}

#endif

// ======================================================================

#include "inc_step.c"

#ifdef PUSHER_BY_BLOCK

static void
do_push_part_1vb_yz(struct psc_fields *flds, struct psc_particles *prts)
{
#ifdef PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);
#endif

  for (int b = 0; b < sub->nr_blocks; b++) {
    for (int n = sub->b_off[b]; n < sub->b_off[b+1]; n++) {
      push_one(flds, prts, n);
    }
  }
}

#else

static void
do_push_part_1vb_yz(struct psc_fields *flds, struct psc_particles *prts)
{
  for (int n = 0; n < prts->n_part; n++) {
    push_one(flds, prts, n);
  }
}

#endif

static void
psc_push_particles_push_a_yz(struct psc_push_particles *push,
			     struct psc_particles *prts,
			     struct psc_fields *flds)
{
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_fields *flds_cache = cache_fields_from_em(flds);

  params_1vb_set(ppsc, flds->p);
  ext_prepare_sort_before(prts);
  do_push_part_1vb_yz(flds_cache, prts);

  cache_fields_to_j(flds_cache, flds);
  psc_fields_destroy(flds_cache);
}

