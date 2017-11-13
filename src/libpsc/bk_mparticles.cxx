
#include "bk_mparticles_iface.h"

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
  return new bk_mparticles;
}

// ----------------------------------------------------------------------
// bk_mparticles_destroy

void bk_mparticles_destroy(struct bk_mparticles *bkmprts)
{
  delete bkmprts;
}

// ----------------------------------------------------------------------
// bk_mparticles_ctor

void bk_mparticles_ctor(struct bk_mparticles *bkmprts, int n_patches)
{
  bkmprts->n_patches = n_patches;

  bkmprts->buf = (PARTICLE_BUF(t) *) calloc(n_patches, sizeof(*bkmprts->buf));
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

void bk_mparticles_reserve_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    PARTICLE_BUF(reserve)(&bkmprts->buf[p], n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_resize_all

void bk_mparticles_resize_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    PARTICLE_BUF(resize)(&bkmprts->buf[p], n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_get_size_all

void bk_mparticles_get_size_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
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



