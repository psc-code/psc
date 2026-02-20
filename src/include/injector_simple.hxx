
#pragma once

#include "particles.hxx"

// ======================================================================
// InjectorSimple

template <typename Mparticles>
struct InjectorSimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;

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

      auto prt = Particle{
        {real_t(new_prt.x[0] - patch.xb[0]), real_t(new_prt.x[1] - patch.xb[1]),
         real_t(new_prt.x[2] - patch.xb[2])},
        {real_t(new_prt.u[0]), real_t(new_prt.u[1]), real_t(new_prt.u[2])},
        real_t(new_prt.w * mprts_.grid().kinds[new_prt.kind].q),
        new_prt.kind,
        mprts_.uid_gen(),
        new_prt.tag};
      mprts_.push_back(p_, prt);
    }

    void reweight(const psc::particle::Inject& new_prt)
    {
      auto& grid = mprts_.grid();
      real_t dVi =
        1.f / (grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2]);
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
