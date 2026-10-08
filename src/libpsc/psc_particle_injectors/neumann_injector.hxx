#pragma once

#include <vector>

#include "grid.hxx"
#include "rng.hxx"
#include "particle.h"
#include <psc.hxx>
#include "pushp.hxx"
#include "kg/Vec3.h"
#include "../psc_bnd/psc_bnd_util.hxx"
#include "injector_base.hxx"

using psc::bnd::LoHi;

/// @brief A resampler for use with @ref NeumannBoundaryInjector that leaves
/// injected copies' velocities unchanged.
class NeumannResamplerNone
{
public:
  using Real3 = psc::particle::Inject::Real3;

  Real3 resample(const psc::particle::Inject& prt, int normal_dim)
  {
    return prt.u;
  }
};

/// @brief A resampler for use with @ref NeumannBoundaryInjector that samples
/// injected copies' transverse velocities from a non-drifting Maxwellian.
class NeumannResamplerMaxwellian
{
public:
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  /// @param temperatures the temperature of each kind, indexed by kind
  NeumannResamplerMaxwellian(const Grid_t::Kinds& kinds,
                             std::vector<Real> temperatures)
  {
    assert(temperatures.size() == kinds.size());
    for (int k = 0; k < kinds.size(); k++) {
      vdfs.emplace_back(0.0, sqrt(temperatures[k] / kinds[k].m));
    }
  }

  Real3 resample(const psc::particle::Inject& prt, int normal_dim)
  {
    auto& vdf = vdfs[prt.kind];
    // FIXME should really sample from Maxwell-Juttner
    // this hack interprets v as u to handle rare case when v>1
    // v<<1 => v~= u anyways
    Real3 u = prt.u;
    for (int d = 0; d < 3; d++) {
      if (d != normal_dim) {
        u[d] = vdf.get();
      }
    }
    return u;
  }

private:
  std::vector<rng::Normal<Real>> vdfs;
};

/// @brief A resampler for use with @ref NeumannBoundaryInjector that rotates
/// injected copies' transverse velocities by a random angle, preserving their
/// magnitudes.
class NeumannResamplerRotate
{
public:
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  Real3 resample(const psc::particle::Inject& prt, int normal_dim)
  {
    int d1 = (normal_dim + 1) % 3;
    int d2 = (normal_dim + 2) % 3;

    Real u_transverse = sqrt(sqr(prt.u[d1]) + sqr(prt.u[d2]));
    Real angle = angle_dist.get();

    Real3 u = prt.u;
    u[d1] = u_transverse * cos(angle);
    u[d2] = u_transverse * sin(angle);
    return u;
  }

private:
  rng::Uniform<Real> angle_dist{0.0, 2.0 * M_PI};
};

/// @brief Injects particles on a given boundary such that the particle
/// distribution satisfies a zero-gradient (von Neumann) boundary condition.
/// Whenever a particle moves from the edge cell inwards to a non-edge cell, a
/// copy of it is injected one cell closer to the boundary, as if there had been
/// an identical particle in the ghost cell. Applies to all particle species.
///
/// Must run after the particle push and before particle boundary exchange, so
/// that the pushed particles are still in their original patches.
///
/// The copies' transverse velocities can optionally be resampled. Their normal
/// momenta are then adjusted to preserve normal velocities, so that the copies
/// still enter from the ghost cell.
/// @tparam LOHI whether to inject at the lower or upper boundary
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`, `AdvanceParticle_t`
/// @tparam RESAMPLER type that defines `resample(prt, normal_dim)`, which takes
/// a copy (as a `psc::particle::Inject`) and the index of the normal dimension,
/// and returns a momentum, of which only the transverse components are used; see @ref NeumannResamplerMaxwellian
template <LoHi LOHI, typename PUSH_PARTICLES,
          typename RESAMPLER = NeumannResamplerNone>
class NeumannBoundaryInjector
  : public InjectorBase<typename PUSH_PARTICLES::Mparticles,
                        typename PUSH_PARTICLES::MfieldsState>
{
  static const int INJECT_DIM_IDX_ = 1;

public:
  using PushParticles = PUSH_PARTICLES;
  using Resampler = RESAMPLER;

  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using Current = typename PushParticles::Current;
  using real_t = typename PushParticles::real_t;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;
  using Real3 = Vec3<real_t>;

  static const bool lo = LOHI == LoHi::Lo;

  NeumannBoundaryInjector(Resampler resampler = {}) : resampler{resampler} {}

  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    const Grid_t& grid = mprts.grid();
    auto injectors_by_patch = mprts.injector();
    auto accessor = mprts.accessor();

    Real3 dxi = grid.domain.dx_inv;
    Current current(grid);
    AdvanceParticle_t advance{grid.dt};

    int edge_idx = lo ? 0 : grid.ldims[INJECT_DIM_IDX_] - 1;
    Real3 ghost_offset = Real3{0, 0, 0}.with_component(
      INJECT_DIM_IDX_, (lo ? -1 : 1) * grid.domain.dx[INJECT_DIM_IDX_]);

    for (int p = 0; p < grid.n_patches(); p++) {
      if (!(lo ? grid.atBoundaryLo(p, INJECT_DIM_IDX_)
               : grid.atBoundaryHi(p, INJECT_DIM_IDX_))) {
        continue;
      }

      // collect copies first, since injecting invalidates the accessor
      std::vector<psc::particle::Inject> copies;
      for (auto prt : accessor[p]) {
        int final_idx = (prt.x() * dxi).fint()[INJECT_DIM_IDX_];

        if (final_idx != edge_idx + (lo ? 1 : -1)) {
          // particles can't move more than one cell at a time (v<c<dx/dt),
          // so only need to check particles in the second-innermost layer
          continue;
        }

        Real3 v = advance.calc_v(prt.u());
        Real3 initial_pos = prt.x();
        advance.push_x(initial_pos, v, -1);

        int initial_idx = (initial_pos * dxi).fint()[INJECT_DIM_IDX_];

        if (initial_idx == edge_idx) {
          using InjectReal3 = psc::particle::Inject::Real3;
          copies.emplace_back(InjectReal3(prt.x() + ghost_offset),
                              InjectReal3(prt.u()), prt.w(), prt.kind(),
                              prt.tag());
          resample(copies.back());
        }
      }

      auto&& injector = injectors_by_patch[p];
      auto flds = mflds[p];
      typename Current::fields_t J(flds);

      for (auto& copy : copies) {
        injector.inject_local(copy);

        Real3 v = advance.calc_v(Real3(copy.u));
        Real3 final_pos = Real3(copy.x);
        Real3 initial_pos = final_pos;
        advance.push_x(initial_pos, v, -1);

        Real3 initial_normalized_pos = initial_pos * dxi;
        Real3 final_normalized_pos = final_pos * dxi;
        Int3 initial_idx = initial_normalized_pos.fint();
        Int3 final_idx = final_normalized_pos.fint();

        current.calc_j(J, initial_normalized_pos, final_normalized_pos,
                       final_idx, initial_idx,
                       grid.kinds[copy.kind].q * real_t(copy.w), v);
      }
    }
  }

private:
  void resample(psc::particle::Inject& copy)
  {
    auto u = copy.u;
    auto v_normal = u[INJECT_DIM_IDX_] / sqrt(1 + u.mag2());

    u = resampler.resample(copy, INJECT_DIM_IDX_);
    u[INJECT_DIM_IDX_] = 0;
    u[INJECT_DIM_IDX_] =
      v_normal * sqrt((1 + u.mag2()) / (1 - v_normal * v_normal));
    copy.u = u;
  }

  Resampler resampler;
};
