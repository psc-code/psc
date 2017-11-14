
#include "bk_mparticles_iface.h"

#include <vector>

#define PARTICLE_BUF(x) psc_particle_single_by_kind_buf_ ## x

struct particle_buf_t {
  typedef particle_single_by_kind_t particle_t;

  particle_buf_t() {
    PARTICLE_BUF(ctor)(&z);
  }

  ~particle_buf_t() {
    PARTICLE_BUF(dtor)(&z);
  }

  void reserve(int n) {
    PARTICLE_BUF(reserve)(&z, n);
  }
  
  void resize(int n) {
    PARTICLE_BUF(resize)(&z, n);
  }
  
  int size() {
    return PARTICLE_BUF(size)(&z);
  }

  particle_t& operator[](int n) {
    return *PARTICLE_BUF(at_ptr)(&z, n);
  }

  psc_particle_single_by_kind_buf_t z;

};

// ======================================================================
// bk_mparticles

struct bk_mparticles {
  typedef particle_single_by_kind_t particle_t;
  //typedef psc_particle_single_by_kind_buf_t particle_buf_t;

  int n_patches;
  std::vector<particle_buf_t> buf;
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
  bkmprts->buf.resize(n_patches);
}

// ----------------------------------------------------------------------
// bk_mparticles_dtor

void bk_mparticles_dtor(struct bk_mparticles *bkmprts)
{
}

// ----------------------------------------------------------------------
// bk_mparticles_reserve_all

void bk_mparticles_reserve_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    bkmprts->buf[p].reserve(n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_resize_all

void bk_mparticles_resize_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    bkmprts->buf[p].resize(n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_get_size_all

void bk_mparticles_get_size_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  for (int p = 0; p < bkmprts->n_patches; p++) {
    n_prts_by_patch[p] = bkmprts->buf[p].size();
  }
}

// ----------------------------------------------------------------------
// bk_mparticles_get_n_prts

int bk_mparticles_get_n_prts(struct bk_mparticles *bkmprts)
{
  int n_prts = 0;
  for (int p = 0; p < bkmprts->n_patches; p++) {
    n_prts += bkmprts->buf[p].size();
  }
  return n_prts;
}

// ----------------------------------------------------------------------
// bk_mparticles_at_ptr

particle_single_by_kind_t *
bk_mparticles_at_ptr(struct bk_mparticles *bkmprts, int p, int n)
{
  return &bkmprts->buf[p][n];
}



