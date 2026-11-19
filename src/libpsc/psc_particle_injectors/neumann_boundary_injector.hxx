#pragma once

#include <vector>

#include "grid.hxx"
#include "rng.hxx"
#include "particle.h"
#include <psc.hxx>
#include "pushp.hxx"
#include "kg/Vec3.h"
#include "../psc_bnd/psc_bnd_util.hxx"
#include "boundary_injector_base.hxx"

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

  /// Resampling the transverse momentum changes gamma, and therefore the
  /// normal velocity if the normal momentum were kept. The current deposited
  /// for the copy assumes it moved from the ghost cell into the domain over the
  /// last step, i.e., it pushes the copy back by its velocity. With a different
  /// normal velocity, that back-push may not cross the boundary, which would
  /// violate charge conservation. So the normal momentum is instead adjusted to
  /// preserve the normal velocity.
  Real3 resample(const psc::particle::Inject& prt, int normal_dim)
  {
    Real v_normal = prt.u[normal_dim] / sqrt(1 + prt.u.mag2());

    auto& vdf = vdfs[prt.kind];
    // FIXME should really sample from Maxwell-Juttner
    // this hack interprets v as u to handle rare case when v>1
    // v<<1 => v~= u anyways
    Real3 u{0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (d != normal_dim) {
        u[d] = vdf.get();
      }
    }

    u[normal_dim] = v_normal * sqrt((1 + u.mag2()) / (1 - sqr(v_normal)));
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
/// Whenever a particle moves inwards across the cell face `copy_depth` cells
/// from the boundary, a copy of it is injected `copy_depth` cells closer to the
/// boundary, as if there had been an identical particle in the ghost cell. With
/// `copy_depth` = 1 (the default) that face is the edge cell's inner face; a
/// deeper face copies from a larger, less inbred population. Applies to all
/// particle species.
///
/// Must run after the particle push and before particle boundary exchange, so
/// that the pushed particles are still in their original patches.
///
/// The copies' transverse velocities can optionally be resampled.
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`, `AdvanceParticle_t`
/// @tparam RESAMPLER type that defines `resample(prt, normal_dim)`, which takes
/// a copy (as a `psc::particle::Inject`) and the index of the normal dimension,
/// and returns the copy's new momentum; see @ref NeumannResamplerMaxwellian
template <typename PUSH_PARTICLES, typename RESAMPLER = NeumannResamplerNone>
class NeumannBoundaryInjector
  : public psc::bnd::particle::BoundaryInjectorBase<PUSH_PARTICLES>
{
  using Base = psc::bnd::particle::BoundaryInjectorBase<PUSH_PARTICLES>;
  using Base::d;

public:
  using PushParticles = PUSH_PARTICLES;
  using Resampler = RESAMPLER;

  using typename Base::Current;
  using typename Base::MfieldsState;
  using typename Base::Mparticles;
  using typename Base::Real3;
  using typename Base::real_t;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;

  /// injects at the lower or upper (`lohi`) boundary normal to `d`
  NeumannBoundaryInjector(psc::Axis d, LoHi lohi, Resampler resampler = {},
                          int copy_depth = 1)
    : Base{d, lohi}, resampler{resampler}, copy_depth{copy_depth}
  {
    assert(copy_depth >= 1);
  }

  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    const Grid_t& grid = mprts.grid();
    this->check_axis(grid);
    auto injectors_by_patch = mprts.injector();
    auto accessor = mprts.accessor();

    Real3 dxi = Real3(grid.domain.dx_inv);
    Current current(grid);
    AdvanceParticle_t advance{real_t(grid.dt)};

    assert(copy_depth < grid.ldims[d]);
    int edge_idx = this->edge_idx(grid);
    int inwards = this->inwards();
    // the copy crosses the boundary face as the particle crosses the copy face
    Real3 ghost_offset = Real3{0, 0, 0}.with_component(
      d, -inwards * copy_depth * grid.domain.dx[d]);

    for (int p = 0; p < grid.n_patches(); p++) {
      if (!this->at_boundary(grid, p)) {
        continue;
      }

      // collect copies first, since injecting invalidates the accessor
      std::vector<psc::particle::Inject> copies;
      for (auto prt : accessor[p]) {
        int final_idx = (prt.x() * dxi).fint()[d];

        if (final_idx != edge_idx + inwards * copy_depth) {
          // particles can't move more than one cell at a time (v<c<dx/dt),
          // so only need to check particles just inside the copy face
          continue;
        }

        Real3 v = advance.calc_v(prt.u());
        Real3 initial_pos = prt.x();
        advance.push_x(initial_pos, v, -1);

        int initial_idx = (initial_pos * dxi).fint()[d];

        if (initial_idx == edge_idx + inwards * (copy_depth - 1)) {
          using InjectReal3 = psc::particle::Inject::Real3;
          copies.emplace_back(InjectReal3(prt.x() + ghost_offset),
                              InjectReal3(prt.u()), prt.w(), prt.kind(),
                              prt.tag());
          auto& copy = copies.back();
          copy.u = resampler.resample(copy, d);
        }
      }

      auto&& injector = injectors_by_patch[p];
      auto flds = mflds[p];
      typename Current::fields_t J(flds);

      for (auto& copy : copies) {
        Real3 v = advance.calc_v(Real3(copy.u));
        Real3 initial_pos = Real3(copy.x);
        advance.push_x(initial_pos, v, -1);
        Real3 initial_normalized_pos = initial_pos * dxi;

        this->inject_with_current(grid, injector, current, J, copy,
                                  initial_normalized_pos,
                                  initial_normalized_pos.fint(), v,
                                  grid.kinds[copy.kind].q * real_t(copy.w));
      }
    }
  }

private:
  Resampler resampler;
  int copy_depth;
};
