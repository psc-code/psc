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
#include "boundary_injector_base.hxx"

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
/// @tparam PARTICLE_GENERATOR a type that defines `get(min_pos, pos_range)` and
/// returns an injectable particle within that range of positions (usually a
/// grid cell); see @ref ParticleGeneratorMaxwellian
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`, etc.
template <typename PARTICLE_GENERATOR, typename PUSH_PARTICLES>
class DirichletBoundaryInjector : public BoundaryInjectorBase<PUSH_PARTICLES>
{
  using Base = BoundaryInjectorBase<PUSH_PARTICLES>;
  using Base::d;

public:
  using ParticleGenerator = PARTICLE_GENERATOR;
  using PushParticles = PUSH_PARTICLES;

  using typename Base::Current;
  using typename Base::MfieldsState;
  using typename Base::Mparticles;
  using typename Base::Real3;
  using typename Base::real_t;
  using InterpolateEM_t = typename PushParticles::InterpolateEM_t;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;

  /// injects at the lower or upper (`lohi`) boundary normal to `d`
  DirichletBoundaryInjector(ParticleGenerator particle_generator, psc::Axis d,
                            LoHi lohi, real_t density = 1.0)
    : Base{d, lohi}, particle_generator_{particle_generator}, density{density}
  {}

  /// Injects particles at the boundary as if there were a population of
  /// particles just beyond the edge. The imaginary particle population is
  /// sampled using the given ParticleGenerator.
  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    const Grid_t& grid = mprts.grid();
    this->check_axis(grid);
    auto injectors_by_patch = mprts.injector();
    bool lo = this->lo();

    Real3 dxi = grid.domain.dx_inv;
    Current current(grid);
    bool preaccelerate = preaccelerate_method != PreaccelerateMethod::None &&
                         preaccelerate_time != 0;

    for (int p = 0; p < grid.n_patches(); p++) {
      if (!this->at_boundary(grid, p)) {
        continue;
      }

      auto&& injector = injectors_by_patch[p];
      auto flds = mflds[p];
      typename Current::fields_t J(flds);
      typename InterpolateEM_t::fields_t EM(flds.storage(), flds.ib());
      InterpolateEM_t ip;
      InterpolateEMBoundary_t ip_boundary;
      AdvanceParticle_t advance{grid.dt};

      Int3 start = Int3{0, 0, 0}.with_component(d, lo ? -1 : grid.ldims[d]);
      Int3 stop = grid.ldims.with_component(d, start[d] + 1);

      for (Int3 initial_idx : VecRange(start, stop)) {
        Real3 cell_corner = Real3(initial_idx) * grid.domain.dx;
        int n_prts_to_try_inject =
          get_n_in_cell(density, grid.norm.prts_per_unit_density, true);

        for (int prt_count = 0; prt_count < n_prts_to_try_inject; prt_count++) {
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
              real_t sample_coord = lo ? 0 : grid.ldims[d];
              sample_fields(
                ip_boundary, EM,
                initial_normalized_pos.with_component(d, sample_coord), e_inner,
                h_inner);
            } else {
              // sample interior E and B, which works for all interpolators
              real_t sample_coord = lo ? 0.5 : grid.ldims[d] - 0.5;
              sample_fields(
                ip, EM, initial_normalized_pos.with_component(d, sample_coord),
                e_inner, h_inner);
            }

            real_t dq = .5f * grid.norm.eta * preaccelerate_time * q / m;
            advance.push_p(prt.u, e_inner, h_inner, dq);
          }

          // push along the normal only, with the relativistic velocity
          Real3 v = advance.calc_v(prt.u);
          prt.x[d] += grid.dt * v[d];

          int final_idx = (prt.x * dxi).fint()[d];
          if (lo ? final_idx <= start[d] : final_idx >= start[d]) {
            // don't inject a particle that fails to enter the patch
            continue;
          }

          this->inject_with_current(grid, injector, current, J, prt,
                                    initial_normalized_pos, initial_idx, v,
                                    q * prt.w);
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
      e[d] = d == 0 ? ip.ex(EM) : d == 1 ? ip.ey(EM) : ip.ez(EM);
    } else if (preaccelerate_method == PreaccelerateMethod::AllE) {
      e = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
    } else if (preaccelerate_method == PreaccelerateMethod::AllEH) {
      e = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
      h = {ip.hx(EM), ip.hy(EM), ip.hz(EM)};
    }
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
