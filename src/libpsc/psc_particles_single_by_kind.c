
#include "psc.h"
#include "psc_particles_as_single_by_kind.h"
#include "psc_particles_inc.h"

// ======================================================================
// psc_mparticles: subclass "single_by_kind"

static struct mrc_obj_method psc_mparticles_single_by_kind_methods[] = {
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup_patch

static void
PFX(setup_patch)(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct PFX(patch) *patch = &sub->patch[p];

  PARTICLE_BUF(ctor)(&patch->buf);

  for (int d = 0; d < 3; d++) {
    patch->b_mx[d] = ppsc->patch[p].ldims[d];
    patch->b_dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_destroy_patch

static void
PFX(destroy_patch)(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct PFX(patch) *patch = &sub->patch[p];

  // need to free structures created in ::patch_setup and ::patch_reserve
  PARTICLE_BUF(dtor)(&patch->buf);
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup

static void
PFX(setup)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  psc_mparticles_setup_super(mprts);
  sub->patch = calloc(mprts->nr_patches, sizeof(*sub->patch));

  for (int p = 0; p < mprts->nr_patches; p++) {
    PFX(setup_patch)(mprts, p);
  }
}

static void
PFX(destroy)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    PFX(destroy_patch)(mprts, p);
  }
  free(sub->patch);
}

static void
PFX(reserve_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    PFX(patch_reserve)(mprts, p, n_prts_by_patch[p]);
  }
}

static void
PFX(resize_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct PFX(patch) *patch = &sub->patch[p];
    PARTICLE_BUF(resize)(&patch->buf, n_prts_by_patch[p]);
  }
}

static void
PFX(get_size_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct PFX(patch) *patch = &sub->patch[p];
    n_prts_by_patch[p] = PARTICLE_BUF(size)(&patch->buf);
  }
}

static unsigned int
PFX(get_nr_particles)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  int n_prts = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct PFX(patch) *patch = &sub->patch[p];
    n_prts += PARTICLE_BUF(size)(&patch->buf);
  }
  return n_prts;
}

// ----------------------------------------------------------------------
// psc_mparticles_ops

struct psc_mparticles_ops PFX(ops) = {
  .name                    = PARTICLE_TYPE,
  .size                    = sizeof(struct psc_mparticles_sub),
  .methods                 = PFX(methods),
  .setup                   = PFX(setup),
  .destroy                 = PFX(destroy),
#if 0
  .write                   = PFX(write),
  .read                    = PFX(read),
#endif
  .reserve_all             = PFX(reserve_all),
  .resize_all              = PFX(resize_all),
  .get_size_all            = PFX(get_size_all),
  .get_nr_particles        = PFX(get_nr_particles),
};

