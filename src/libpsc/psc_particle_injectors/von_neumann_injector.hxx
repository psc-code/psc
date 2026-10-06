#pragma once

#include <vector>

#include "grid.hxx"
#include "particle.h"
#include <psc.hxx>
#include "pushp.hxx"
#include "kg/Vec3.h"
#include "../psc_bnd/psc_bnd_util.hxx"
#include "injector_base.hxx"

using psc::bnd::LoHi;

/// @brief Injects particles on a given boundary such that the particle
/// distribution satisfies a zero-gradient (von Neumann) boundary condition.
/// Whenever a particle moves from the edge cell inwards to a non-edge cell, a
/// copy of it is injected one cell closer to the boundary, as if there had been
/// an identical particle in the ghost cell. Applies to all particle species.
///
/// Must run after the particle push and before particle boundary exchange, so
/// that the pushed particles are still in their original patches.
/// @tparam LOHI whether to inject at the lower or upper boundary
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`, `AdvanceParticle_t`
template <LoHi LOHI, typename PUSH_PARTICLES>
class VonNeumannBoundaryInjector
  : public InjectorBase<typename PUSH_PARTICLES::Mparticles,
                        typename PUSH_PARTICLES::MfieldsState>
{
  static const int INJECT_DIM_IDX_ = 1;

public:
  using PushParticles = PUSH_PARTICLES;

  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using Current = typename PushParticles::Current;
  using real_t = typename PushParticles::real_t;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;
  using Real3 = Vec3<real_t>;

  static const bool lo = LOHI == LoHi::Lo;

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
};
