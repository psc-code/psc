
#include "heating.hxx"

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
    float ran1, ran2, ran3, ran4, ran5, ran6;
    do {
      ran1 = random() / ((float)RAND_MAX + 1);
      ran2 = random() / ((float)RAND_MAX + 1);
      ran3 = random() / ((float)RAND_MAX + 1);
      ran4 = random() / ((float)RAND_MAX + 1);
      ran5 = random() / ((float)RAND_MAX + 1);
      ran6 = random() / ((float)RAND_MAX + 1);
    } while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f || ran4 >= 1.f ||
             ran5 >= 1.f || ran6 >= 1.f);

    real_t ranx = sqrtf(-2.f * logf(1.0 - ran1)) * cosf(2.f * M_PI * ran2);
    real_t rany = sqrtf(-2.f * logf(1.0 - ran3)) * cosf(2.f * M_PI * ran4);
    real_t ranz = sqrtf(-2.f * logf(1.0 - ran5)) * cosf(2.f * M_PI * ran6);

    real_t Dpxi = sqrtf(H * heating_dt_);
    real_t Dpyi = sqrtf(H * heating_dt_);
    real_t Dpzi = sqrtf(H * heating_dt_);

    prt.u[0] += Dpxi * ranx;
    prt.u[1] += Dpyi * rany;
    prt.u[2] += Dpzi * ranz;
  }

  void operator()(Mparticles& mprts)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto&& prts = mprts[p];
      auto& patch = mprts.grid().patches[p];
      for (auto prt_iter = mprts.begin(p); prt_iter != mprts.end(p);
           prt_iter++) {
        double xx[3] = {
          prt_iter->x[0] + patch.xb[0],
          prt_iter->x[1] + patch.xb[1],
          prt_iter->x[2] + patch.xb[2],
        };
        double H = get_H_(xx, prt_iter->kind);
        if (H > 0.f) {
          kick_particle(*prt_iter, H);
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
