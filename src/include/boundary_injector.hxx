#pragma once

#include <math.h>

#include "grid.hxx"
#include "rng.hxx"
#include "particle.h"
#include <psc.hxx>
#include "pushp.hxx"
#include "dim.hxx"
#include "setup_particles.hxx"
#include "../libpsc/psc_push_particles/inc_push.cxx"

class ParticleGeneratorMaxwellian
{
public:
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  // FIXME would be nice to just pass 1 thing for kind-related info
  ParticleGeneratorMaxwellian(int kind_idx, Grid_t::Kind kind, Real3 mean_u,
                              Real3 temperature)
    : kind_idx{kind_idx}
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
    Real w = 1.0;
    psc::particle::Tag tag = 0;

    return {x, u, w, kind_idx, tag};
  }

private:
  using VelocityDistributionFunction = rng::Normal<Real>;
  Vec3<VelocityDistributionFunction> vdfs;
  int kind_idx;
  rng::Uniform<Real> uniform_dist{0.0, 1.0};
};

template <typename PARTICLE_GENERATOR, typename PUSH_PARTICLES>
class BoundaryInjector
{
  static const int INJECT_DIM_IDX_ = 1;

public:
  using ParticleGenerator = PARTICLE_GENERATOR;
  using PushParticles = PUSH_PARTICLES;

  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;
  using Current = typename PushParticles::Current;
  using Dim = typename PushParticles::Dim;
  using real_t = typename PushParticles::real_t;
  using Real3 = typename PushParticles::Real3;
  using checks_order = typename PushParticles::checks_order;

  BoundaryInjector(ParticleGenerator particle_generator, Grid_t& grid)
    : partice_generator{particle_generator},
      advance{grid.dt},
      prts_per_unit_density{grid.norm.prts_per_unit_density}
  {}

  /// Injects particles at the lower y-bound as if there were a population of
  /// particles just beyond the edge. The imaginary particle population has unit
  /// density, and individual particles from that population are sampled using
  /// the given ParticleGenerator.
  ///
  /// Some of these limitations may be removed in the future.
  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    static_assert(INJECT_DIM_IDX_ == 1,
                  "only injection at lower bound of y is supported");

    const Grid_t& grid = mprts.grid();
    auto injectors_by_patch = mprts.injector();

    Real3 dxi = Real3{1., 1., 1.} / Real3(grid.domain.dx);
    Current current(grid);

    for (int patch_idx = 0; patch_idx < grid.n_patches(); patch_idx++) {
      Int3 ilo = grid.patches[patch_idx].off;

      if (ilo[INJECT_DIM_IDX_] != 0) {
        continue;
      }

      Int3 ihi = ilo + grid.ldims;

      auto&& injector = injectors_by_patch[patch_idx];
      auto flds = mflds[patch_idx];
      typename Current::fields_t J(flds);

      for (Int3 cell_idx = ilo; cell_idx[0] < ihi[0]; cell_idx[0]++) {
        for (cell_idx[2] = ilo[2]; cell_idx[2] < ihi[2]; cell_idx[2]++) {
          auto boundary_current_before =
            flds.storage()(cell_idx[0] + grid.ibn[0], -1 + grid.ibn[1],
                           cell_idx[2] + grid.ibn[2], JXI + 1);

          cell_idx[INJECT_DIM_IDX_] = -1;

          Real3 cell_corner =
            grid.domain.corner + Double3(cell_idx) * grid.domain.dx;
          int n_prts_to_try_inject =
            get_n_in_cell(1.0, prts_per_unit_density, true);
          real_t charge_injected = 0;

          for (int prt_count = 0; prt_count < n_prts_to_try_inject;
               prt_count++) {
            psc::particle::Inject prt =
              partice_generator.get(cell_corner, grid.domain.dx);

            Real3 v = advance.calc_v(prt.u);
            Real3 initial_x = prt.x;
            advance.push_x(prt.x, v, 1.0);

            if (prt.x[INJECT_DIM_IDX_] < grid.domain.corner[INJECT_DIM_IDX_]) {
              // don't inject a particle that fails to enter the domain
              continue;
            }

            injector(prt);

            // Update currents
            // Taken from push_particles_1vb.hxx PushParticlesVb::push_mprts()

            Real3 initial_normalized_pos = initial_x * dxi;

            Int3 final_idx;
            Real3 final_normalized_pos;
            find_idx_pos_1st_rel(prt.x, dxi, final_idx, final_normalized_pos,
                                 real_t(0.));

            // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
            real_t qni_wni = grid.kinds[prt.kind].q * prt.w;
            current.calc_j(J, initial_normalized_pos, final_normalized_pos,
                           final_idx, cell_idx, qni_wni, v);

            charge_injected += qni_wni;
          }

          // override whatever current was deposited in the boundary to be
          // consistent with the charge having started completely out of the
          // domain
          flds.storage()(cell_idx[0] + grid.ibn[0], -1 + grid.ibn[1],
                         cell_idx[2] + grid.ibn[2], JXI + 1) =
            boundary_current_before + charge_injected * grid.domain.dx[1] /
                                        grid.dt / prts_per_unit_density;
        }
      }
    }
  }

private:
  ParticleGenerator partice_generator;
  AdvanceParticle<real_t, dim_y> advance;
  real_t prts_per_unit_density;
};
