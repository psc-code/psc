
#include "heating.hxx"
#include "rng.hxx"

#include <functional>
#include <stdlib.h>

// ======================================================================
// Heating__

template <typename MP>
struct Heating__ : HeatingBase
{
  using Mparticles = MP;
  using real_t = typename Mparticles::real_t;
  using Particle = typename Mparticles::Particle;

  // ----------------------------------------------------------------------
  // ctor

  template <typename FUNC>
  Heating__(const Grid_t& grid, int interval, FUNC get_H) : get_H_{get_H}
  {
    heating_dt_ = interval * grid.dt;
  }

  // ----------------------------------------------------------------------
  // kick_particle

  void kick_particle(Particle& prt, real_t H)
  {
    rng::Normal<real_t> standard_normal_distr;

    real_t Dpxi = sqrtf(H * heating_dt_);
    real_t Dpyi = sqrtf(H * heating_dt_);
    real_t Dpzi = sqrtf(H * heating_dt_);

    prt.u[0] += Dpxi * standard_normal_distr.get();
    prt.u[1] += Dpyi * standard_normal_distr.get();
    prt.u[2] += Dpzi * standard_normal_distr.get();
  }

  void operator()(Mparticles& mprts)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto&& prts = mprts[p];
      auto& patch = mprts.grid().patches[p];
      for (auto& prt : prts) {

        double xx[3] = {
          prt.x[0] + patch.xb[0],
          prt.x[1] + patch.xb[1],
          prt.x[2] + patch.xb[2],
        };
        double H = get_H_(xx, prt.kind);
        if (H > 0.f) {
          kick_particle(prt, H);
        }
      }
    }
  }

private:
  real_t heating_dt_;
  std::function<double(const double*, const int)> get_H_;
};

// ======================================================================
// HeatingSelector
//
// FIXME, this should become unnecessary

template <typename Mparticles>
struct HeatingSelector
{
  using Heating = Heating__<Mparticles>;
};

#ifdef USE_CUDA

#include "../libpsc/cuda/heating_cuda_impl.hxx"

// FIXME, enable_if for any BS
template <>
struct HeatingSelector<MparticlesCuda<BS444>>
{
  using Mparticles = MparticlesCuda<BS444>;
  using Heating = HeatingCuda<HeatingSpotFoil<dim_xyz>, Mparticles>;
};

template <>
struct HeatingSelector<MparticlesCuda<BS144>>
{
  using Mparticles = MparticlesCuda<BS144>;
  using Heating = HeatingCuda<HeatingSpotFoil<dim_yz>, Mparticles>;
};

#endif
