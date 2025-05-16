#pragma once

#include <math.h>

#include "grid.hxx"
#include "rng.hxx"
#include "particle.h"
#include "pushp.hxx"
#include "dim.hxx"

class ParticleGeneratorMaxwellian
{
public:
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  // FIXME would be nice to just pass 1 thing for kind-related info
  ParticleGeneratorMaxwellian(int kind_idx, Grid_t::Kind kind, Real w,
                              Real3 mean_u, Real3 temperature)
    : kind_idx{kind_idx}, w{w}
  {
    for (int d = 0; d < 3; d++) {
      Real stdev_u = sqrt(temperature[d] / kind.m);
      vdfs[d] = VelocityDistributionFunction{mean_u[d], stdev_u};
    }
  }

  psc::particle::Inject get(Real3 min_pos, Real3 pos_range)
  {
    Real3 x;
    for (int d = 0; d < 3; d++) {
      x[d] = min_pos[d] + uniform_dist.get() * pos_range[d];
    }

    Real3 u{vdfs[0].get(), vdfs[1].get(), vdfs[2].get()};
    psc::particle::Tag tag = 0;

    return {x, u, w, kind_idx, tag};
  }

private:
  using VelocityDistributionFunction = rng::Normal<Real>;
  Vec3<VelocityDistributionFunction> vdfs;
  Real w;
  int kind_idx;
  rng::Uniform<Real> uniform_dist{0.0, 1.0};
};

template <typename PARTICLE_GENERATOR, typename PUSH_PARTICLES>
class InjectorBoundaryInflow
{
  static const int INJECT_DIM_IDX_ = 1;

public:
  using ParticleGenerator = PARTICLE_GENERATOR;
  using PushParticles = PUSH_PARTICLES;

  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;
  using InterpolateEM_t = typename PushParticles::InterpolateEM_t;
  using Current = typename PushParticles::Current;
  using Dim = typename PushParticles::Dim;
  using real_t = typename PushParticles::real_t;
  using Real3 = typename PushParticles::Real3;
  using checks_order = typename PushParticles::checks_order;

  InjectorBoundaryInflow(ParticleGenerator particle_generator, Grid_t& grid)
    : partice_generator{particle_generator},
      advance{grid.dt},
      prts_per_cell{grid.norm.prts_per_unit_density}
  {}

  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    const Grid_t& grid = mprts.grid();
    auto injectors_by_patch = mprts.injector();

    for (int patch_idx = 0; patch_idx < grid.n_patches(); patch_idx++) {
      const auto& patch = grid.patches[patch_idx];

      if (patch.off[INJECT_DIM_IDX_] != 0) {
        continue;
      }

      auto&& injector = injectors_by_patch[patch_idx];

      Int3 ilo = patch.off;
      Int3 ihi = ilo + grid.ldims;

      assert(INJECT_DIM_IDX_ == 1);
      assert(ilo[INJECT_DIM_IDX_] == 0);

      for (Int3 cell_idx = ilo; cell_idx[0] < ihi[0]; cell_idx[0]++) {
        for (cell_idx[2] = ilo[2]; cell_idx[2] < ihi[2]; cell_idx[2]++) {
          assert(cell_idx[INJECT_DIM_IDX_] == 0);
          cell_idx[INJECT_DIM_IDX_] = -1;

          for (int prt_count = 0; prt_count < prts_per_cell; prt_count++) {
            Real3 min_pos =
              grid.domain.corner + Double3(cell_idx) * grid.domain.dx;

            psc::particle::Inject prt =
              partice_generator.get(min_pos, grid.domain.dx);

            Real3 v = advance.calc_v(prt.u);
            advance.push_x(prt.x, v, 1.0);

            if (prt.x[INJECT_DIM_IDX_] < grid.domain.corner[INJECT_DIM_IDX_]) {
              continue;
            }

            injector(prt);
          }
        }
      }
    }

    // TODO:
    // 4. update current
  }

private:
  ParticleGenerator partice_generator;
  AdvanceParticle<real_t, dim_y> advance;
  real_t prts_per_cell;
};
