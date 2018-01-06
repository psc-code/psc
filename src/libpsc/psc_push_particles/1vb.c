
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define USE_1VB

#include "fields.hxx"

#include "inc_params.c"
#include "inc_cache.c"
#include "inc_interpolate.c"
#include "inc_push.c"
#include "inc_curr.c"
#include "inc_step.c"

// ======================================================================

#ifdef PUSHER_BY_BLOCK

static void
do_push_part_1vb_yz(fields_t flds, struct psc_mparticles *mprts, int p)
{
#ifdef PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  struct psc_mparticles_single_by_block *msub = psc_mparticles_single_by_block(mprts);
  struct psc_mparticles_single_by_block_patch *patch = &msub->patch[p];
#endif

  particle_range_t prts = particle_range_mprts(mprts, p);
  for (int b = 0; b < patch->nr_blocks; b++) {
    for (int n = patch->b_off[b]; n < patch->b_off[b+1]; n++) {
      push_one(prts.begin, n, flds, flds);
    }
  }
}

static void
do_stagger_part_1vb_yz(fields_t flds, struct psc_mparticles *mprts, int p)
{
#ifdef PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  struct psc_mparticles_single_by_block *msub = psc_mparticles_single_by_block(mprts);
  struct psc_mparticles_single_by_block_patch *patch = &msub->patch[p];
#endif

  particle_range_t prts = particle_range_mprts(mprts, p);
  for (int b = 0; b < patch->nr_blocks; b++) {
    for (int n = patch->b_off[b]; n < patch->b_off[b+1]; n++) {
      stagger_one(prts.begin, n, flds);
    }
  }
}

#else

static void
do_push_part_1vb_yz(fields_t flds, struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = particle_range_mprts(mprts, p);
  unsigned int n_prts = particle_range_size(prts);
  
  for (int n = 0; n < n_prts; n++) {
    push_one(prts.begin, n, flds, flds);
  }
}

static void
do_stagger_part_1vb_yz(fields_t flds, struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = particle_range_mprts(mprts, p);
  unsigned int n_prts = particle_range_size(prts);
  
  for (int n = 0; n < n_prts; n++) {
    stagger_one(prts.begin, n, flds);
  }
}

#endif

void
SFX(psc_push_particles_push_mprts)(struct psc_push_particles *push,
				   struct psc_mparticles *mprts,
				   struct psc_mfields *mflds)
{
  mfields_t mf(mflds);
  c_prm_set(ppsc);
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = mf[p];

    fields_t_zero_range(flds, JXI, JXI + 3);
    ext_prepare_sort_before(mprts, p);
    do_push_part_1vb_yz(flds, mprts, p);
  }
}

void
SFX(psc_push_particles_stagger_mprts)(struct psc_push_particles *push,
				      struct psc_mparticles *mprts,
				      struct psc_mfields *mflds)
{
  mfields_t mf(mflds);
  c_prm_set(ppsc);
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = mf[p];

    fields_t_zero_range(flds, JXI, JXI + 3);
    ext_prepare_sort_before(mprts, p);
    do_stagger_part_1vb_yz(flds, mprts, p);
  }
}

