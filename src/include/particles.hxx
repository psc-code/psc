
#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include "psc_particles.h"

// ======================================================================
// mparticles_base

template<typename S>
struct mparticles_base
{
  using sub_t = S;

  explicit mparticles_base(psc_mparticles *mprts) : mprts_(mprts) { }

  void put_as(psc_mparticles *mprts_base, unsigned int flags = 0)
  {
    psc_mparticles_put_as(mprts_, mprts_base, flags);
  }

  unsigned int n_patches() { return mprts_->nr_patches; }
  
  psc_mparticles *mprts() { return mprts_; }
  
  sub_t* sub() { return mrc_to_subobj(mprts(), sub_t); }

private:
  psc_mparticles *mprts_;
};

// ======================================================================
// psc_particle

template<class R>
struct psc_particle
{
  using real_t = R;

  real_t xi, yi, zi;
  real_t qni_wni;
  real_t pxi, pyi, pzi;
  int kind_;

  int kind() { return kind_; }
};

// ======================================================================
// mparticles

template<typename S>
struct mparticles : mparticles_base<S>
{
  using Base = mparticles_base<S>;
  using particles_t = typename Base::sub_t::particles_t;
  using particle_buf_t = typename particles_t::buf_t;
  
  mparticles(psc_mparticles *mprts) : mparticles_base<S>(mprts) { }

  particles_t& operator[](int p)
  {
    return this->sub()->patch[p];
  }

};

#endif

