
#include "bk_mparticles_iface.h"

// ----------------------------------------------------------------------
// mparticles::ctor

template<typename T>
mparticles<T>::mparticles(int _n_patches)
  : n_patches(_n_patches), buf(_n_patches)
{
}

// ----------------------------------------------------------------------
// mparticles::reserve_all

template<typename T>
void mparticles<T>::reserve_all(const int n_prts_by_patch[])
{
  for (int p = 0; p < n_patches; p++) {
    buf[p].reserve(n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// mparticles::resize_all

template<typename T>
void mparticles<T>::resize_all(const int n_prts_by_patch[])
{
  for (int p = 0; p < n_patches; p++) {
    buf[p].resize(n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// mparticles::size_all

template<typename T>
void mparticles<T>::size_all(int n_prts_by_patch[]) const
{
  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = buf[p].size();
  }
}

// ----------------------------------------------------------------------
// mparticles::n_prts

template<typename T>
int mparticles<T>::n_prts() const
{
  int n_prts = 0;
  for (int p = 0; p < n_patches; p++) {
    n_prts += buf[p].size();
  }
  return n_prts;
}

// ----------------------------------------------------------------------
// mparticles::at

template<typename T>
typename mparticles<T>::particle_type& mparticles<T>::at(int p, int n)
{
  return buf[p][n];
}

// ======================================================================
// bk_mparticles C wrappers

bk_mparticles *bk_mparticles_new(int n_patches)
{
  return new bk_mparticles(n_patches);
}

void bk_mparticles_delete(bk_mparticles *bkmprts)
{
  delete bkmprts;
}

void bk_mparticles_reserve_all(bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  bkmprts->reserve_all(n_prts_by_patch);
}

void bk_mparticles_resize_all(bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  bkmprts->resize_all(n_prts_by_patch);
}

void bk_mparticles_size_all(bk_mparticles *bkmprts, int n_prts_by_patch[])
{
  bkmprts->size_all(n_prts_by_patch);
}

int bk_mparticles_n_prts(bk_mparticles *bkmprts)
{
  return bkmprts->n_prts();
}

particle_single_by_kind_t *bk_mparticles_at_ptr(bk_mparticles *bkmprts, int p, int n)
{
  return &bkmprts->at(p, n);
}



