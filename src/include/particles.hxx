
#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include "psc_particles.h"
#include "grid.hxx"

#include <vector>

// ======================================================================
// psc_particle

template<class R>
struct psc_particle
{
  using real_t = R;

  real_t xi, yi, zi;
  real_t qni_wni_;
  real_t pxi, pyi, pzi;
  int kind_;

  int kind() { return kind_; }

  // FIXME, grid is always double precision, so this will switch precision
  // where not desired. should use same info stored in mprts at right precision
  real_t qni(const Grid_t& grid) const { return grid.kinds[kind_].q; }
  real_t mni(const Grid_t& grid) const { return grid.kinds[kind_].m; }
  real_t wni(const Grid_t& grid) const { return qni_wni_ / qni(grid); }
  real_t qni_wni(const Grid_t& grid) const { return qni_wni_; }
};

// ======================================================================
// mparticles_base

template<typename S>
struct mparticles_base
{
  using sub_t = S;
  using particle_t = typename sub_t::particle_t;
  using real_t = typename particle_t::real_t;

  explicit mparticles_base(psc_mparticles *mprts)
    : mprts_(mprts),
      sub_(mrc_to_subobj(mprts, sub_t))
  {
  }

  void put_as(psc_mparticles *mprts_base, unsigned int flags = 0)
  {
    psc_mparticles_put_as(mprts_, mprts_base, flags);
    mprts_ = nullptr;
  }

  void get_size_all(int *n_prts_by_patch)
  {
    psc_mparticles_get_size_all(mprts_, n_prts_by_patch);
  }

  void reserve_all(int *n_prts_by_patch)
  {
    psc_mparticles_reserve_all(mprts_, n_prts_by_patch);
  }

  void resize_all(int *n_prts_by_patch)
  {
    psc_mparticles_resize_all(mprts_, n_prts_by_patch);
  }

  unsigned int n_patches() { return mprts_->nr_patches; }
  
  psc_mparticles *mprts() { return mprts_; }
  
  sub_t* operator->() { return sub_; }

public:
  psc_mparticles *mprts_;
  sub_t *sub_;
};

// ======================================================================
// mparticles_patch_base

template<typename P>
struct psc_mparticles_;

template<typename P>
struct mparticles_patch_base
{
  using particle_t = P;
  using real_t = typename particle_t::real_t;
  using buf_t = std::vector<particle_t>;
  using iterator = typename buf_t::iterator;
  
  buf_t buf = {};

  int b_mx[3];
  real_t b_dxi[3];

  psc_mparticles_<P>* mprts;
  int p;

  // FIXME, I would like to delete the copy ctor because I don't
  // want to copy patch_t by mistake, but that doesn't play well with
  // putting the patches into std::vector
  // mparticles_patch_base(const mparticles_patch_base&) = delete;

  mparticles_patch_base(psc_mparticles_<P>* mprts, int p)
  {
    this->mprts = mprts;
    this->p = p;
    
    for (int d = 0; d < 3; d++) {
      b_mx[d] = mprts->grid_.ldims[d];
      b_dxi[d] = 1.f / mprts->grid_.dx[d];
    }
    
#if PSC_PARTICLES_AS_SINGLE
    patch->nr_blocks = patch->b_mx[0] * patch->b_mx[1] * patch->b_mx[2];
    patch->b_cnt = (unsigned int *) calloc(patch->nr_blocks + 1, sizeof(*patch->b_cnt));
#endif
  }

  particle_t& operator[](int n) { return buf[n]; }
  iterator begin() { return buf.begin(); }
  iterator end() { return buf.end(); }
  unsigned int size() const { return buf.size(); }
  void reserve(unsigned int new_capacity) { buf.reserve(new_capacity); }
  void push_back(const particle_t& prt) { buf.push_back(prt); }

  void resize(unsigned int new_size)
  {
    assert(new_size <= buf.capacity());
    buf.resize(new_size);
  }

  buf_t& get_buf()
  {
    return buf;
  }

  void get_block_pos(const real_t xi[3], int b_pos[3])
  {
    for (int d = 0; d < 3; d++) {
      b_pos[d] = fint(xi[d] * b_dxi[d]);
    }
  }
  
  const int* get_b_mx() const { return b_mx; }
  const real_t* get_b_dxi() const { return b_dxi; }
};

template<typename P>
struct mparticles_patch : mparticles_patch_base<P> {
  using Base = mparticles_patch_base<P>;
  
  using Base::Base;
};

// ======================================================================
// psc_mparticles_

template<typename P>
struct psc_mparticles_
{
  using particle_t = P;
  using particle_real_t = typename particle_t::real_t; // FIXME, should go away
  using real_t = particle_real_t;
  using patch_t = mparticles_patch<particle_t>;

  psc_mparticles_(const Grid_t& grid)
    : grid_(grid)
  {
    patches_.reserve(grid.patches.size());
    for (int p = 0; p < grid.patches.size(); p++) {
      patches_.emplace_back(this, p);
    }
  }

  const patch_t& operator[](int p) const { return patches_[p]; }
  patch_t&       operator[](int p)       { return patches_[p]; }
  
  particle_real_t prt_qni(const particle_t& prt) const { return prt.qni(grid_); }
  particle_real_t prt_mni(const particle_t& prt) const { return prt.mni(grid_); }
  particle_real_t prt_wni(const particle_t& prt) const { return prt.wni(grid_); }
  particle_real_t prt_qni_wni(const particle_t& prt) const { return prt.qni_wni(grid_); }

  std::vector<patch_t> patches_;
  const Grid_t& grid_;
};

// ======================================================================
// mparticles

template<typename S>
struct mparticles : mparticles_base<S>
{
  using Base = mparticles_base<S>;
  using patch_t = typename Base::sub_t::patch_t;
  using particle_buf_t = typename patch_t::buf_t;
  
  mparticles(psc_mparticles *mprts) : mparticles_base<S>(mprts) { }

  patch_t& operator[](int p)
  {
    return (*this->sub_)[p];
  }

};

#endif

