
#pragma once

#include "psc_particles.h"

// ======================================================================
// InjectorBuffered

template<typename Mparticles>
struct InjectorBuffered
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = typename Particle::Real3;
  using Double3 = Vec3<double>;
  
  struct Patch
  {
    Patch(InjectorBuffered& injector, int p)
      : injector_(injector), p_{p}, n_prts_{0} // FIXME, why (), {} does not work
    {}
    
    ~Patch()
    {
      injector_.n_prts_by_patch_[p_] += n_prts_;
    }
  
    void raw(const Particle& prt)
    {
      injector_.buf_.push_back(prt);
      n_prts_++;
    }
    
    void operator()(const psc::particle::Inject& new_prt)
    {
      auto& mprts = injector_.mprts_;
      auto& patch = mprts.grid().patches[p_];
      auto x = Double3::fromPointer(new_prt.x) - patch.xb;
      auto u = Double3::fromPointer(new_prt.u);
      real_t q = mprts.grid().kinds[new_prt.kind].q;
      injector_.buf_.push_back({Real3(x), Real3(u), q * real_t(new_prt.w), new_prt.kind, mprts.uid_gen(), new_prt.tag});
      n_prts_++;
    }

    // FIXME do we want to keep this? or just have a psc::particle::Inject version instead?
    void raw(const std::vector<Particle>& buf)
    {
      injector_.buf_.insert(injector_.buf_.end(), buf.begin(), buf.end());
      n_prts_ += buf.size();
    }
    
  private:
    InjectorBuffered& injector_;
    const int p_;
    uint n_prts_;
  };
  
  InjectorBuffered(Mparticles& mprts)
    : n_prts_by_patch_(mprts.n_patches()), last_patch_{0}, mprts_{mprts}
  {}

  ~InjectorBuffered()
  {
    assert(n_prts_by_patch_.size() == mprts_.n_patches());
    mprts_.inject(buf_, n_prts_by_patch_);
  }

  Patch operator[](int p)
  {
    // ensure that we inject particles into patches in ascending order
    assert(p >= last_patch_);
    last_patch_ = p;
    return {*this, p};
  }

private:
  std::vector<Particle> buf_;
  std::vector<uint> n_prts_by_patch_;
  uint last_patch_;
  Mparticles& mprts_;
};

