#pragma once

#include <math.h>

#include "grid.hxx"
#include "rng.hxx"
#include "particle.h"
#include <psc.hxx>
#include "pushp.hxx"
#include "dim.hxx"
#include "setup_particles.hxx"
#include "kg/VecRange.hxx"
#include "../libpsc/psc_push_particles/inc_push.cxx"
#include "injector_base.hxx"

/// @brief A particle generator for use with @ref BoundaryInjector. Samples
/// particles from a (possibly shifted) Maxwellian distribution.
class ParticleGeneratorMaxwellian
{
public:
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  // FIXME would be nice to just pass 1 thing for kind-related info
  ParticleGeneratorMaxwellian(int kind_idx, Grid_t::Kind kind, Real3 mean_v,
                              Real3 temperature)
    : kind_idx{kind_idx}, prt_booster{-mean_v}
  {
    for (int d = 0; d < 3; d++) {
      vdfs[d] = VDF{0.0, sqrt(temperature[d] / kind.m)};
    }
  }

  psc::particle::Inject get(Real3 min_pos, Real3 pos_range)
  {
    Real3 x;
    for (int d = 0; d < 3; d++) {
      x[d] = min_pos[d] + uniform_dist.get() * pos_range[d];
    }

    Real3 v{vdfs[0].get(), vdfs[1].get(), vdfs[2].get()};
    // FIXME should really sample from Maxwell-Juttner
    // this hack interprests v as u to handle rare case when v>1
    // v<<1 => v~= u anyways
    Real3 u = prt_booster.boost(v);

    Real w = 1.0;
    psc::particle::Tag tag = 0;

    return {x, u, w, kind_idx, tag};
  }

private:
  using VDF = rng::Normal<Real>;
  Vec3<VDF> vdfs;
  VelocityBooster prt_booster;
  int kind_idx;
  rng::Uniform<Real> uniform_dist{0.0, 1.0};
};

/// @brief Injects particles on a given boundary, sampling from a given particle
/// generator. For precise control over multiple particle species, use one
/// BoundaryInjector per species.
/// @tparam PARTICLE_GENERATOR a type that defines `get(min_pos, pos_range)` and
/// returns an injectable particle within that range of positions (usually a
/// grid cell); see @ref ParticleGeneratorMaxwellian
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`, etc.
template <typename PARTICLE_GENERATOR, typename PUSH_PARTICLES>
class BoundaryInjector
  : public InjectorBase<typename PUSH_PARTICLES::Mparticles,
                        typename PUSH_PARTICLES::MfieldsState>
{
  static const int INJECT_DIM_IDX_ = 1;

public:
  using ParticleGenerator = PARTICLE_GENERATOR;
  using PushParticles = PUSH_PARTICLES;

  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using Current = typename PushParticles::Current;
  using real_t = typename PushParticles::real_t;
  using InterpolateEM_t = typename PushParticles::InterpolateEM_t;
  using Real3 = Vec3<real_t>;

  BoundaryInjector(ParticleGenerator particle_generator)
    : particle_generator_{particle_generator}
  {}

  /// Injects particles at specified y-bounds as if there were a population of
  /// particles just beyond the edge. The imaginary particle population is
  /// sampled using the given ParticleGenerator.
  ///
  /// The dimensional limitations may be removed in the future.
  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    static_assert(INJECT_DIM_IDX_ == 1,
                  "only injection at lower bound of y is supported");

    const Grid_t& grid = mprts.grid();
    auto injectors_by_patch = mprts.injector();

    Real3 dxi = grid.domain.dx_inv;
    Current current(grid);

    for (int p = 0; p < grid.n_patches(); p++) {
      // TODO: combine paths as much as possible
      if (inject_lo && grid.atBoundaryLo(p, INJECT_DIM_IDX_)) {
        Int3 ilo = {0, 0, 0};
        Int3 ihi = grid.ldims;

        ilo[INJECT_DIM_IDX_] = -1;
        ihi[INJECT_DIM_IDX_] = 0;

        auto&& injector = injectors_by_patch[p];
        auto flds = mflds[p];
        typename Current::fields_t J(flds);
        typename InterpolateEM_t::fields_t EM(flds.storage(), flds.ib());
        InterpolateEM_t ip;
        AdvanceParticle<real_t, dim_y> advance{grid.dt};

        for (Int3 initial_idx : VecRange(ilo, ihi)) {
          Real3 cell_corner = Real3(initial_idx) * grid.domain.dx;
          int n_prts_to_try_inject =
            get_n_in_cell(density, grid.norm.prts_per_unit_density, true);

          for (int prt_count = 0; prt_count < n_prts_to_try_inject;
               prt_count++) {
            // 1. sample position uniformly from first ghost layer, and velocity
            // from vdf at x=infty
            psc::particle::Inject prt =
              particle_generator_.get(cell_corner, grid.domain.dx);

            real_t m = grid.kinds[prt.kind].m;
            real_t q = grid.kinds[prt.kind].q;

            Real3 initial_normalized_pos = prt.x * dxi;

            bool preaccelerate = true;
            if (preaccelerate) {
              // 2. calculate change in energy from traveling from -infty to xp
              real_t xp = initial_normalized_pos[INJECT_DIM_IDX_];
              // x0 < xp < x1 and E0=E(x0), E1=E(x1)
              real_t x0, x1;
              real_t E0, E1;

              if (xp < -0.5) {
                x0 = -1.5;
                E0 = 0.0;

                x1 = -0.5;
              } else {
                x0 = -0.5;
                ip.set_coeffs(
                  initial_normalized_pos.with_component(INJECT_DIM_IDX_, x0));
                switch (INJECT_DIM_IDX_) {
                  case 0: E0 = ip.ex(EM); break;
                  case 1: E0 = ip.ey(EM); break;
                  case 2: E0 = ip.ez(EM); break;
                  default: assert(false);
                }

                x1 = 0.5;
              }
              ip.set_coeffs(
                initial_normalized_pos.with_component(INJECT_DIM_IDX_, x1));
              switch (INJECT_DIM_IDX_) {
                case 0: E1 = ip.ex(EM); break;
                case 1: E1 = ip.ey(EM); break;
                case 2: E1 = ip.ez(EM); break;
                default: assert(false);
              }

              real_t dx = grid.domain.dx[INJECT_DIM_IDX_];
              real_t work =
                0.5f * q * dx *
                (E0 * (2 - sqr(xp - x0 - 1)) + E1 * sqr(xp - x1 + 1));

              // 3. recalculate normal velocity based on work already done
              real_t u_inf = prt.u[INJECT_DIM_IDX_];
              real_t gamma_inf = sqrt(1 + prt.u.mag2());
              real_t gamma_inf_perp =
                sqrt(1 + prt.u.with_component(INJECT_DIM_IDX_, 0).mag2());
              real_t gamma_p = gamma_inf + work / m;
              if (gamma_p < gamma_inf_perp) {
                // particle didn't have enough energy to reach xp
                continue;
              }
              // sign is a hack
              prt.u[INJECT_DIM_IDX_] =
                sqrt(sqr(gamma_p) - sqr(gamma_inf_perp)) * (u_inf < 0 ? -1 : 1);
              prt.w *= u_inf * gamma_p / (prt.u[INJECT_DIM_IDX_] * gamma_inf);

              // 4. push normal u a half step
              real_t E_interp = E0 * (1 - xp + x0) + E1 * (1 + xp - x1);
              prt.u[INJECT_DIM_IDX_] += 0.5f * grid.dt * q * E_interp / m;
            }

            // 5. push normal x
            Real3 v = advance.calc_v(prt.u);
            advance.push_x(prt.x, v);

            if (prt.x[INJECT_DIM_IDX_] < 0.0) {
              // don't inject a particle that fails to enter the patch
              continue;
            }

            Real3 final_normalized_pos = prt.x * dxi;
            Int3 final_idx = final_normalized_pos.fint();

            injector.inject_local(prt);

            real_t qni_wni = grid.kinds[prt.kind].q * prt.w;
            current.calc_j(J, initial_normalized_pos, final_normalized_pos,
                           final_idx, initial_idx, qni_wni, v);
          }
        }
      }

      if (inject_hi && grid.atBoundaryHi(p, INJECT_DIM_IDX_)) {
        Int3 ilo = {0, 0, 0};
        Int3 ihi = grid.ldims;

        ilo[INJECT_DIM_IDX_] = grid.ldims[INJECT_DIM_IDX_];
        ihi[INJECT_DIM_IDX_] = grid.ldims[INJECT_DIM_IDX_] + 1;

        auto&& injector = injectors_by_patch[p];
        auto flds = mflds[p];
        typename Current::fields_t J(flds);

        for (Int3 initial_idx : VecRange(ilo, ihi)) {
          Real3 cell_corner = Real3(initial_idx) * grid.domain.dx;
          int n_prts_to_try_inject =
            get_n_in_cell(density, grid.norm.prts_per_unit_density, true);

          Int3 inner_idx = initial_idx - Int3::unit(INJECT_DIM_IDX_);
          real_t inner_e = mflds(EX + INJECT_DIM_IDX_, inner_idx[0],
                                 inner_idx[1], inner_idx[2], p);

          for (int prt_count = 0; prt_count < n_prts_to_try_inject;
               prt_count++) {
            psc::particle::Inject prt =
              particle_generator_.get(cell_corner, grid.domain.dx);

            // FIXME no factor of dt - why does this work? can it be better?
            real_t du_from_e =
              inner_e * grid.kinds[prt.kind].q / grid.kinds[prt.kind].m;
            prt.u[INJECT_DIM_IDX_] += du_from_e;

            AdvanceParticle<real_t, dim_y> advance{grid.dt};
            Real3 v = advance.calc_v(prt.u);
            Real3 initial_normalized_pos = prt.x * dxi;
            advance.push_x(prt.x, v);
            Real3 final_normalized_pos = prt.x * dxi;
            Int3 final_idx = final_normalized_pos.fint();

            if (final_idx[INJECT_DIM_IDX_] >= initial_idx[INJECT_DIM_IDX_]) {
              // don't inject a particle that fails to enter the patch
              continue;
            }

            injector.inject_local(prt);

            real_t qni_wni = grid.kinds[prt.kind].q * prt.w;
            current.calc_j(J, initial_normalized_pos, final_normalized_pos,
                           final_idx, initial_idx, qni_wni, v);
          }
        }
      }
    }
  }

public:
  real_t density = 1.0;
  bool inject_lo = true;
  bool inject_hi = false;

private:
  ParticleGenerator particle_generator_;
};
