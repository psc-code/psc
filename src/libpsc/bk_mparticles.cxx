
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
  
  int size() const {
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

  void reserve_all(const int n_prts_by_patch[]);
  void resize_all(const int n_prts_by_patch[]);
  void size_all(int n_prts_by_patch[]) const;
  int n_prts() const;
  particle_t& at(int p, int n);
};

// ----------------------------------------------------------------------
// bk_mparticles::reserve_all

void bk_mparticles::reserve_all(const int n_prts_by_patch[])
{
  for (int p = 0; p < n_patches; p++) {
    buf[p].reserve(n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles::resize_all

void bk_mparticles::resize_all(const int n_prts_by_patch[])
{
  for (int p = 0; p < n_patches; p++) {
    buf[p].resize(n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// bk_mparticles::size_all

void bk_mparticles::size_all(int n_prts_by_patch[]) const
{
  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = buf[p].size();
  }
}

// ----------------------------------------------------------------------
// bk_mparticles::n_prts

int bk_mparticles::n_prts() const
{
  int n_prts = 0;
  for (int p = 0; p < n_patches; p++) {
    n_prts += buf[p].size();
  }
  return n_prts;
}

// ----------------------------------------------------------------------
// bk_mparticles::at

bk_mparticles::particle_t& bk_mparticles::at(int p, int n)
{
  return buf[p][n];
}

// ======================================================================
// bk_mparticles C wrappers

struct bk_mparticles *bk_mparticles_create()
{
  return new bk_mparticles;
}

void bk_mparticles_destroy(struct bk_mparticles *bkmprts)
{
  delete bkmprts;
}

void bk_mparticles_ctor(struct bk_mparticles *bkmprts, int n_patches)
{
  bkmprts->n_patches = n_patches;
  bkmprts->buf.resize(n_patches);
}

void bk_mparticles_dtor(struct bk_mparticles *bkmprts)
{
}

void bk_mparticles_reserve_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  bkmprts->reserve_all(n_prts_by_patch);
}

void bk_mparticles_resize_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  bkmprts->resize_all(n_prts_by_patch);
}

void bk_mparticles_size_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  bkmprts->size_all(n_prts_by_patch);
}

int bk_mparticles_n_prts(struct bk_mparticles *bkmprts)
{
  return bkmprts->n_prts();
}

particle_single_by_kind_t *bk_mparticles_at_ptr(struct bk_mparticles *bkmprts, int p, int n)
{
  return &bkmprts->at(p, n);
}



