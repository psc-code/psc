
#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include "psc_particles.h"

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

template<typename P, typename S>
struct mparticles : mparticles_base<S>
{
  using Base = mparticles_base<S>;
  using mparticles_patch_t = typename Base::sub_t::patch_t;
  using buf_t = P;
  using particle_t = typename buf_t::particle_t;
  
  mparticles(psc_mparticles *mprts) : mparticles_base<S>(mprts) { }

  mparticles_patch_t& operator[](int p)
  {
    return this->sub()->patch[p];
  }

};

#endif

