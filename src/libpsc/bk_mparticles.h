
#ifndef BK_MPARTICLES_H
#define BK_MPARTICLES_H

#include <vector>

// ======================================================================
// bk_mparticles

template<typename T>
struct mparticles {
  typedef T particle_type;
  typedef std::vector<particle_type> particle_buf_t;

  int n_patches;
  std::vector<particle_buf_t> buf;

  mparticles(int _n_patches);

  void reserve_all(const int n_prts_by_patch[]);
  void resize_all(const int n_prts_by_patch[]);
  void size_all(int n_prts_by_patch[]) const;
  int n_prts() const;
  particle_type& at(int p, int n);
};

typedef mparticles<particle_single_by_kind_t> bk_mparticles;

#endif
