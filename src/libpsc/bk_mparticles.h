
#ifndef BK_MPARTICLES_H
#define BK_MPARTICLES_H

#include <vector>

// ======================================================================
// mparticles

template<typename T>
struct mparticles {
  typedef T particle_buf_type;
  typedef typename T::value_type particle_type;

  int n_patches;
  std::vector<particle_buf_type> buf;

  mparticles(int _n_patches);

  void reserve_all(const int n_prts_by_patch[]);
  void resize_all(const int n_prts_by_patch[]);
  void size_all(int n_prts_by_patch[]) const;
  int n_prts() const;
  particle_type& at(int p, int n);
};

#endif
