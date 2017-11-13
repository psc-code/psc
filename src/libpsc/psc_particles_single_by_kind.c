
#include "psc.h"
#include "psc_particles_single_by_kind.h"

#define PARTICLE_TYPE "single_by_kind"
#define PFX(x) psc_mparticles_single_by_kind_ ## x
#define psc_mparticles_sub psc_mparticles_single_by_kind
#define PARTICLE_BUF(x) psc_particle_single_by_kind_buf_ ## x

// ======================================================================
// bk_mparticles

struct bk_mparticles {
  int n_patches;
  PARTICLE_BUF(t) **buf;
};

// ----------------------------------------------------------------------
// bk_mparticles_create

struct bk_mparticles *bk_mparticles_create()
{
  return calloc(1, sizeof(struct bk_mparticles));
}

// ----------------------------------------------------------------------
// bk_mparticles_ctor

void bk_mparticles_ctor(struct bk_mparticles *bkmprts, int n_patches,
			PARTICLE_BUF(t) **_buf)
{
  bkmprts->n_patches = n_patches;
  bkmprts->buf = calloc(n_patches, sizeof(*bkmprts->buf));
}

// ----------------------------------------------------------------------
// bk_mparticles_dtor

void bk_mparticles_dtor(struct bk_mparticles *bkmprts)
{
  free(bkmprts->buf);
}

// ======================================================================
// psc_mparticles: subclass "single_by_kind"

struct psc_mparticles_single_by_kind_patch {
  psc_particle_single_by_kind_buf_t buf;

  int b_mx[3];
  particle_single_by_kind_real_t b_dxi[3];
};

static struct mrc_obj_method psc_mparticles_single_by_kind_methods[] = {
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup

static void
PFX(setup)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  psc_mparticles_setup_super(mprts);
  sub->patch = calloc(mprts->nr_patches, sizeof(*sub->patch));

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct PFX(patch) *patch = &sub->patch[p];
    
    PARTICLE_BUF(ctor)(&patch->buf);
    
    for (int d = 0; d < 3; d++) {
      patch->b_mx[d] = ppsc->patch[p].ldims[d];
      patch->b_dxi[d] = 1.f / ppsc->patch[p].dx[d];
    }
  }

  sub->bkmprts = bk_mparticles_create(); // FIXME, leaked

  PARTICLE_BUF(t) *buf[mprts->nr_patches];
  for (int p = 0; p < mprts->nr_patches; p++) {
    buf[p] = &sub->patch[p].buf;
  }
  bk_mparticles_ctor(sub->bkmprts, mprts->nr_patches, buf);
}

static void
PFX(destroy)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct PFX(patch) *patch = &sub->patch[p];
    
    // need to free structures created in ::patch_setup and ::patch_reserve
    PARTICLE_BUF(dtor)(&patch->buf);
  }
  free(sub->patch);

  bk_mparticles_dtor(sub->bkmprts);
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

void
psc_mparticles_single_by_kind_patch_reserve(struct psc_mparticles *mprts, int p,
					    unsigned int new_capacity)
{
  struct psc_mparticles_single_by_kind *sub = psc_mparticles_single_by_kind(mprts);
  struct psc_mparticles_single_by_kind_patch *patch = &sub->patch[p];

  psc_particle_single_by_kind_buf_reserve(&patch->buf, new_capacity);
}

particle_single_by_kind_t *
psc_mparticles_single_by_kind_get_one(struct psc_mparticles *mprts, int p, unsigned int n)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_single_by_kind_ops);
  struct psc_mparticles_single_by_kind *sub = psc_mparticles_single_by_kind(mprts);
  struct psc_mparticles_single_by_kind_patch *patch = &sub->patch[p];

  return psc_particle_single_by_kind_buf_at_ptr(&patch->buf, n);
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

