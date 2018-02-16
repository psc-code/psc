
#ifndef BK_MPARTICLES_H
#define BK_MPARTICLES_H

#include <psc_bits.h>
#include <vector>

// ======================================================================
// mparticles

template<typename T>
struct mparticles_ {
  typedef T particle_buf_type;
  typedef typename T::value_type particle_type;

  int n_patches;
  std::vector<particle_buf_type> buf;

  // ----------------------------------------------------------------------
  // mparticles::ctor

  mparticles_(int _n_patches)
    : n_patches(_n_patches), buf(_n_patches)
  {
  }

  // ----------------------------------------------------------------------
  // mparticles::reserve_all

  void reserve_all(const uint n_prts_by_patch[])
  {
    for (int p = 0; p < n_patches; p++) {
      buf[p].reserve(n_prts_by_patch[p]);
    }
  }

  // ----------------------------------------------------------------------
  // mparticles::resize_all
  
  void resize_all(const uint n_prts_by_patch[])
  {
    for (int p = 0; p < n_patches; p++) {
      buf[p].resize(n_prts_by_patch[p]);
    }
  }

  // ----------------------------------------------------------------------
  // mparticles::size_all
  
  void size_all(uint n_prts_by_patch[]) const
  {
    for (int p = 0; p < n_patches; p++) {
      n_prts_by_patch[p] = buf[p].size();
    }
  }

  // ----------------------------------------------------------------------
  // mparticles::n_prts
  
  int n_prts() const
  {
    int n_prts = 0;
    for (int p = 0; p < n_patches; p++) {
      n_prts += buf[p].size();
    }
    return n_prts;
  }

  // ----------------------------------------------------------------------
  // mparticles::at
  
  particle_type& at(int p, int n)
  {
    return buf[p][n];
  }
};

#endif
