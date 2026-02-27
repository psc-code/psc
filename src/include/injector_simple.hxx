
#pragma once

#include "particles.hxx"
#include "kg/Vec3.h"

// ======================================================================
// InjectorSimple

template <typename Mparticles>
struct InjectorSimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = Vec3<real_t>;

  struct Patch
  {
    Patch(Mparticles& mprts, int p) : mprts_{mprts}, p_{p} {}

    void operator()(const psc::particle::Inject& new_prt)
    {
      const auto& patch = mprts_.grid().patches[p_];
      for (int d = 0; d < 3; d++) {
        assert(new_prt.x[d] >= patch.xb[d]);
        assert(new_prt.x[d] <= patch.xe[d]);
      }

      auto prt =
        Particle{Real3(new_prt.x) - Real3(patch.xb),
                 Real3(new_prt.u),
                 real_t(new_prt.w * mprts_.grid().kinds[new_prt.kind].q),
                 new_prt.kind,
                 mprts_.uid_gen(),
                 new_prt.tag};
      mprts_.push_back(p_, prt);
    }

    void reweight(const psc::particle::Inject& new_prt)
    {
      auto& grid = mprts_.grid();
      real_t dVi = 1.f / grid.domain.dx.prod();
      auto prt = new_prt;
      assert(0); // FIXME, have to actually do reweighting
      (*this)(prt);
    }

  private:
    Mparticles& mprts_;
    int p_;
  };

  InjectorSimple(Mparticles& mprts) : mprts_{mprts} {}

  void reserve(int n_prts_total) {}

  Patch operator[](int p) const { return {mprts_, p}; }

private:
  Mparticles& mprts_;
};
