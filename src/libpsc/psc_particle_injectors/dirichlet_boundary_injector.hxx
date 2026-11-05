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
#include "../psc_push_particles/inc_push.cxx"
#include "../psc_bnd/psc_bnd_util.hxx"
#include "injector_base.hxx"

using psc::bnd::LoHi;

enum class PreaccelerateMethod
{
  None,
  NormalE,
  AllE,
  AllEH,
};

/// @brief A particle generator for use with @ref DirichletBoundaryInjector.
/// Samples particles from a (possibly shifted) Maxwellian distribution.
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
/// DirichletBoundaryInjector per species.
/// @tparam LOHI whether to inject at the lower or upper boundary
/// @tparam PARTICLE_GENERATOR a type that defines `get(min_pos, pos_range)` and
/// returns an injectable particle within that range of positions (usually a
/// grid cell); see @ref ParticleGeneratorMaxwellian
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`, etc.
template <LoHi LOHI, typename PARTICLE_GENERATOR, typename PUSH_PARTICLES>
class DirichletBoundaryInjector
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

  static const bool lo = LOHI == LoHi::Lo;

  DirichletBoundaryInjector(ParticleGenerator particle_generator,
                            real_t density = 1.0)
    : particle_generator_{particle_generator}, density{density}
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
    bool preaccelerate = preaccelerate_method != PreaccelerateMethod::None &&
                         preaccelerate_time != 0;

    for (int p = 0; p < grid.n_patches(); p++) {
      if (lo ? grid.atBoundaryLo(p, INJECT_DIM_IDX_)
             : grid.atBoundaryHi(p, INJECT_DIM_IDX_)) {
        auto&& injector = injectors_by_patch[p];
        auto flds = mflds[p];
        typename Current::fields_t J(flds);
        typename InterpolateEM_t::fields_t EM(flds.storage(), flds.ib());
        InterpolateEM_t ip;
        InterpolateEMBoundary_t ip_boundary;
        AdvanceParticle<real_t, dim_y> advance{grid.dt};

        Int3 start = Int3{0, 0, 0}.with_component(
          INJECT_DIM_IDX_, lo ? -1 : grid.ldims[INJECT_DIM_IDX_]);
        Int3 stop = grid.ldims.with_component(INJECT_DIM_IDX_,
                                              start[INJECT_DIM_IDX_] + 1);

        for (Int3 initial_idx : VecRange(start, stop)) {
          Real3 cell_corner = Real3(initial_idx) * grid.domain.dx;
          int n_prts_to_try_inject =
            get_n_in_cell(density, grid.norm.prts_per_unit_density, true);

          for (int prt_count = 0; prt_count < n_prts_to_try_inject;
               prt_count++) {
            // sample position uniformly from first ghost layer, and velocity
            // from vdf at x=infty
            psc::particle::Inject prt =
              particle_generator_.get(cell_corner, grid.domain.dx);

            real_t m = grid.kinds[prt.kind].m;
            real_t q = grid.kinds[prt.kind].q;

            Real3 initial_normalized_pos = prt.x * dxi;

            if (preaccelerate) {
              Real3 e_inner = {0.f, 0.f, 0.f};
              Real3 h_inner = {0.f, 0.f, 0.f};

              if (preaccelerate_sample_boundary) {
                // the standard 1st-order interpolator really interpolates
                // across the boundary, averaging ghost and edge cells
                real_t sample_coord = lo ? 0 : grid.ldims[INJECT_DIM_IDX_];
                sample_fields(ip_boundary, EM,
                              initial_normalized_pos.with_component(
                                INJECT_DIM_IDX_, sample_coord),
                              e_inner, h_inner);
              } else {
                // sample interior E and B, which works for all interpolators
                real_t sample_coord =
                  lo ? 0.5 : grid.ldims[INJECT_DIM_IDX_] - 0.5;
                sample_fields(ip, EM,
                              initial_normalized_pos.with_component(
                                INJECT_DIM_IDX_, sample_coord),
                              e_inner, h_inner);
              }

              real_t dq = .5f * grid.norm.eta * preaccelerate_time * q / m;
              advance.push_p(prt.u, e_inner, h_inner, dq);
            }

            // push normal x
            Real3 v = advance.calc_v(prt.u);
            advance.push_x(prt.x, v);

            Real3 final_normalized_pos = prt.x * dxi;
            Int3 final_idx = final_normalized_pos.fint();

            if (lo ? final_idx[INJECT_DIM_IDX_] <= start[INJECT_DIM_IDX_]
                   : final_idx[INJECT_DIM_IDX_] >= start[INJECT_DIM_IDX_]) {
              // don't inject a particle that fails to enter the patch
              continue;
            }

            injector.inject_local(prt);

            // deposit as if the particle came from the outer edge of the
            // ghost cell, so its virtual charge lies entirely beyond the
            // boundary node and Gauss's law holds there
            Real3 virtual_initial_pos = initial_normalized_pos.with_component(
              INJECT_DIM_IDX_, outer_ghost_edge(grid));
            current.calc_j(J, virtual_initial_pos, final_normalized_pos,
                           final_idx, initial_idx, q * prt.w, v);
          }
        }
      }
    }
  }

private:
  using InterpolateEMBoundary_t =
    InterpolateEM<typename InterpolateEM_t::fields_t, opt_ip_1st,
                  typename PushParticles::Dim>;

  /// fields at normalized position `pos`, as selected by preaccelerate_method
  template <typename IP, typename F>
  void sample_fields(IP& ip, const F& EM, Real3 pos, Real3& e, Real3& h) const
  {
    ip.set_coeffs(pos);
    if (preaccelerate_method == PreaccelerateMethod::NormalE) {
      e[INJECT_DIM_IDX_] = INJECT_DIM_IDX_ == 0   ? ip.ex(EM)
                           : INJECT_DIM_IDX_ == 1 ? ip.ey(EM)
                                                  : ip.ez(EM);
    } else if (preaccelerate_method == PreaccelerateMethod::AllE) {
      e = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
    } else if (preaccelerate_method == PreaccelerateMethod::AllEH) {
      e = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
      h = {ip.hx(EM), ip.hy(EM), ip.hz(EM)};
    }
  }

  /// normalized coordinate of the ghost cell's face away from the domain
  real_t outer_ghost_edge(const Grid_t& grid) const
  {
    return lo ? -1 : grid.ldims[INJECT_DIM_IDX_] + 1;
  }

public:
  real_t density;
  real_t preaccelerate_time = 0.0;
  PreaccelerateMethod preaccelerate_method = PreaccelerateMethod::NormalE;
  /// sample the preaccelerating fields at the boundary itself, averaging the
  /// ghost and edge cells, instead of at the center of the edge cell
  bool preaccelerate_sample_boundary = false;

private:
  ParticleGenerator particle_generator_;
};
