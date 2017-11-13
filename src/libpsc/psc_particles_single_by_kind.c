
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
  psc_particle_single_by_kind_buf_t *buf;
};

// ----------------------------------------------------------------------
// bk_mparticles_create

struct bk_mparticles *bk_mparticles_create()
{
  return calloc(1, sizeof(struct bk_mparticles));
}

// ----------------------------------------------------------------------
// bk_mparticles_ctor

void bk_mparticles_ctor(struct bk_mparticles *bkmprts, int n_patches)
{
  bkmprts->n_patches = n_patches;

  bkmprts->buf = calloc(n_patches, sizeof(*bkmprts->buf));
  for (int p = 0; p < n_patches; p++) {
    PARTICLE_BUF(ctor)(&bkmprts->buf[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_dtor

void bk_mparticles_dtor(struct bk_mparticles *bkmprts)
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    // need to free structures created in ::patch_setup and ::patch_reserve
    PARTICLE_BUF(dtor)(&bkmprts->buf[p]);
  }
  free(bkmprts->buf);
}

// ----------------------------------------------------------------------
// bk_mparticles_reserve_all

void bk_mparticles_reserve_all(struct bk_mparticles *bkmprts,
			       int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    PARTICLE_BUF(reserve)(&bkmprts->buf[p], n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_resize_all

void bk_mparticles_resize_all(struct bk_mparticles *bkmprts,
			       int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    PARTICLE_BUF(resize)(&bkmprts->buf[p], n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_get_size_all

void bk_mparticles_get_size_all(struct bk_mparticles *bkmprts,
				int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    n_prts_by_patch[p] = PARTICLE_BUF(size)(&bkmprts->buf[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_get_n_prts

int bk_mparticles_get_n_prts(struct bk_mparticles *bkmprts)
{
  int n_prts = 0;
  for (int p = 0; p < bkmprts->n_patches; p++) {
    n_prts += PARTICLE_BUF(size)(&bkmprts->buf[p]);
  }
  return n_prts;
}

// ----------------------------------------------------------------------
// bk_mparticles_at_ptr

particle_single_by_kind_t *
bk_mparticles_at_ptr(struct bk_mparticles *bkmprts, int p, int n)
{
  return PARTICLE_BUF(at_ptr)(&bkmprts->buf[p], n);
}

// ======================================================================
// psc_mparticles: subclass "single_by_kind"

static struct mrc_obj_method psc_mparticles_single_by_kind_methods[] = {
  {}
};

void psc_mparticles_single_by_kind_patch_reserve(struct psc_mparticles *mprts, int p,
						 unsigned int new_capacity);

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup

static void
PFX(setup)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  psc_mparticles_setup_super(mprts);

  sub->bkmprts = bk_mparticles_create();
  bk_mparticles_ctor(sub->bkmprts, mprts->nr_patches);
}

static void
PFX(destroy)(struct psc_mparticles *mprts)
{
  struct bk_mparticles *bkmprts = psc_mparticles_sub(mprts)->bkmprts;

  bk_mparticles_dtor(bkmprts);
}

static void
PFX(reserve_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct bk_mparticles *bkmprts = psc_mparticles_sub(mprts)->bkmprts;

  bk_mparticles_reserve_all(bkmprts, n_prts_by_patch);
}

static void
PFX(resize_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct bk_mparticles *bkmprts = psc_mparticles_sub(mprts)->bkmprts;

  bk_mparticles_resize_all(bkmprts, n_prts_by_patch);
}

static void
PFX(get_size_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct bk_mparticles *bkmprts = psc_mparticles_sub(mprts)->bkmprts;

  bk_mparticles_get_size_all(bkmprts, n_prts_by_patch);
}

static unsigned int
PFX(get_nr_particles)(struct psc_mparticles *mprts)
{
  struct bk_mparticles *bkmprts = psc_mparticles_sub(mprts)->bkmprts;

  return bk_mparticles_get_n_prts(bkmprts);
}

particle_single_by_kind_t *
PFX(get_one)(struct psc_mparticles *mprts, int p, unsigned int n)
{
  struct bk_mparticles *bkmprts = psc_mparticles_sub(mprts)->bkmprts;

  return bk_mparticles_at_ptr(bkmprts, p, n);
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

