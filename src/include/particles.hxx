
#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include "psc_particles.h"
#include "grid.hxx"
#include "psc_bits.h"

#include <vector>


// ======================================================================
// Particle positions in cell / patch

// Particle positions are stored as patch-relative positions, however,
// it is required to know the exact cell a particle is in in a number of places:
// - for performance sorting particles
// - for keeping particles sorted in, e.g., the CUDA particle pusher
// - for finding the appropriate EM fields to interpolate
// - for correctly depositing currents
// - for collisions
//
// The goal here is to establish rules / invariants of the position of
// particles to where (in which patch) they are stored and how to
// recover the cell they are in.
//
// To complicate things, there are currently two aspects to this: Cell
// position and block position, where the former refers to the
// computational mesh that E, B and J live on, whereas a block refers
// to a fixed-size super-cell (e.g., 4x4x4 cells), motivated by
// performance considerations. It is currently not necessarily clear
// that the calculated block indices and cell indices satisfy the
// anticipated relation (bpos[d] = cpos[d] / bs[d]) because of potential
// finite precision arithmetic
//
// Rules / invariants:
//
// for all particles in a given patch,
// (1) the cell position cpos
//     calculated in a given dimension d will satisfy 0 <= cpos < ldims[d]
// (2) the patch relative position xi satisfies
//     0 <= xi <= xm[d] = ldims[d] * dx[d]
//     with the equality at the upper limit only being allowed at a right/top
//     non-periodic boundary
//
// These invariants will be temporarily violated after the particle push, but will
// be restored by the bnd exchange.
//
// Tricky issues to deal with:
// - (1) and (2) should be closely related, but finite precision
//   arithmetic can cause surprises.
//
//   E.g.: dx = 1 ldims = 100. Particle ends up at position -1e-7. It
//   gets sent to the left, where it's new position will be -1e-7 +
//   100., which is actually = 100. (in single precision), meaning that
//   particle violates (1) and (2) in its new patch.
//
// - Calculating cpos correctly when legally xi == ldims[d] * xi[d] at a right boundary
//
// TODO:
// - have cell index be the primary quantity computed, always, and find
//   block index from that
// - boundary exchange should be based on cell, not block index
//   (though if the two indices are always consistent, it becomes a non-issue)

// ======================================================================
// ParticleIndexer

template<class R>
struct ParticleIndexer
{
  using real_t = R;
  using Real3 = Vec3<real_t>;

  ParticleIndexer(const Grid_t& grid)
    : dxi_(Real3(1.) / Real3(grid.dx)),
      ldims_(grid.ldims),
      b_dxi_(Real3(1.) / (Real3(grid.bs) * Real3(grid.dx))),
      b_mx_(grid.ldims / grid.bs)
  {
    n_cells_ = ldims_[0] * ldims_[1] * ldims_[2];
  }

  int cellPosition(real_t x, int d) const
  {
    return fint(x * dxi_[d]);
  }

  int blockPosition(real_t x, int d) const
  {
    return fint(x * b_dxi_[d]);
  }

  Int3 cellPosition(const real_t* pos) const
  {
    Int3 idx;
    for (int d = 0; d < 3; d++) {
      idx[d] = cellPosition(pos[d], d);
    }
    return idx;
  }

  Int3 blockPosition(const real_t* pos) const
  {
    Int3 idx;
    for (int d = 0; d < 3; d++) {
      idx[d] = blockPosition(pos[d], d);
    }
    return idx;
  }

  int cellIndex(const Int3& cpos) const
  {
    if (uint(cpos[0]) >= ldims_[0] ||
	uint(cpos[1]) >= ldims_[1] ||
	uint(cpos[2]) >= ldims_[2]) {
      return -1;
    }
    
    return (cpos[2] * ldims_[1] + cpos[1]) * ldims_[0] + cpos[0];
  }
  
  int blockIndex(const Int3& bpos) const
  {
    if (uint(bpos[0]) >= b_mx_[0] ||
	uint(bpos[1]) >= b_mx_[1] ||
	uint(bpos[2]) >= b_mx_[2]) {
      return -1;
    } else {
      return (bpos[2] * b_mx_[1] + bpos[1]) * b_mx_[0] + bpos[0];
    }
  }
  
  int cellIndex(const real_t* pos) const
  {
    Int3 cpos = cellPosition(pos);
    return cellIndex(cpos);
  }

  int validCellIndex(const real_t* pos) const
  {
    Int3 cpos = cellPosition(pos);
    for (int d = 0; d < 3; d++) {
      if (uint(cpos[d]) >= ldims_[d]) {
	printf("validCellIndex: cpos[%d] = %d ldims_[%d] = %d // pos[%d] = %g pos[%d]*dxi_[%d] = %g\n",
	       d, cpos[d], d, ldims_[d], d, pos[d], d, d, pos[d] * dxi_[d]);
	assert(0);
      }
    }
    int cidx = cellIndex(cpos);
    assert(cidx >= 0);
    return cidx;
  }

  int blockIndex(const real_t* pos) const
  {
    Int3 bpos = blockPosition(pos);
    return blockIndex(bpos);
  }

  void checkInPatchMod(real_t* pos) const
  {
    for (int d = 0; d < 3; d++) {
      int bpos = blockPosition(pos[d], d);
      if (bpos < 0 || bpos >= b_mx_[d]) {
	printf("checkInPatchMod pos %g %g %g\n", pos[0], pos[1], pos[2]);
	printf("checkInPatchMod d %d pos %g bpos %d // %d\n",
	       d, pos[d], bpos, b_mx_[d]);
	if (bpos < 0) {
	  pos[d] = 0.f;
	} else {
	  pos[d] *= (1. - 1e-6);
	}
	bpos = blockPosition(pos[d], d);
      }
      assert(bpos >= 0 && bpos < b_mx_[d]);
    }
  } 
  
  //private:
  Real3 dxi_;
  Int3 ldims_;
  uint n_cells_;
  Real3 b_dxi_;
  Int3 b_mx_;
};

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
  using Real3 = Vec3<real_t>;
  using buf_t = std::vector<particle_t>;
  using iterator = typename buf_t::iterator;
  
  buf_t buf;

  ParticleIndexer<real_t> pi_;

  psc_mparticles_<P>* mprts;
  int p;

  // FIXME, I would like to delete the copy ctor because I don't
  // want to copy patch_t by mistake, but that doesn't play well with
  // putting the patches into std::vector
  // mparticles_patch_base(const mparticles_patch_base&) = delete;

  mparticles_patch_base(psc_mparticles_<P>* _mprts, int _p)
    : pi_(_mprts->grid_),
      mprts(_mprts),
      p(_p)
  {}

  particle_t& operator[](int n) { return buf[n]; }
  iterator begin() { return buf.begin(); }
  iterator end() { return buf.end(); }
  unsigned int size() const { return buf.size(); }
  void reserve(unsigned int new_capacity) { buf.reserve(new_capacity); }

  void push_back(particle_t& prt) // FIXME, should particle_t be const?
  {
    checkInPatchMod(prt);
    validCellIndex(prt);
    buf.push_back(prt);
  }

  void resize(unsigned int new_size)
  {
    assert(new_size <= buf.capacity());
    buf.resize(new_size);
  }

  buf_t& get_buf()
  {
    return buf;
  }

  void check() const
  {
    for (auto& prt : buf) {
      validCellIndex(prt);
    }
  }

  // ParticleIndexer functionality
  int cellPosition(real_t xi, int d) const { return pi_.cellPosition(xi, d); }
  int validCellIndex(const particle_t& prt) const { return pi_.validCellIndex(&prt.xi); }

  int blockPosition(real_t xi, int d) const { return pi_.blockPosition(xi, d); }
  Int3 blockPosition(const real_t* xi) const { return pi_.blockPosition(xi); }
  int blockIndex(const real_t* xi) const { return pi_.blockIndex(xi); }
  void checkInPatchMod(particle_t& prt) const { return pi_.checkInPatchMod(&prt.xi); }
    
  const int* get_b_mx() const { return pi_.b_mx_; }
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
    patches_.reserve(grid.n_patches());
    for (int p = 0; p < grid.n_patches(); p++) {
      patches_.emplace_back(this, p);
    }
  }

  const patch_t& operator[](int p) const { return patches_[p]; }
  patch_t&       operator[](int p)       { return patches_[p]; }

  void reserve_all(const int *n_prts_by_patch)
  {
    for (int p = 0; p < patches_.size(); p++) {
      patches_[p].reserve(n_prts_by_patch[p]);
    }
  }

  void resize_all(const int *n_prts_by_patch)
  {
    for (int p = 0; p < patches_.size(); p++) {
      patches_[p].resize(n_prts_by_patch[p]);
    }
  }

  void get_size_all(int *n_prts_by_patch) const
  {
    for (int p = 0; p < patches_.size(); p++) {
      n_prts_by_patch[p] = patches_[p].size();
    }
  }

  int get_n_prts() const
  {
    int n_prts = 0;
    for (auto const& patch : patches_) {
      n_prts += patch.size();
    }
    return n_prts;
  }
  
  void check() const
  {
    for (auto& patch: patches_) {
      patch.check();
    }
  }
  
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

