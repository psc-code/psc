#pragma once

#include "grid.hxx"
#include "injector_base.hxx"

/// @brief Completes the current of particles leaving through open boundaries.
///
/// After the push, an outgoing particle sits in the first ghost cell, at a
/// depth set by its speed, until particle boundary exchange deletes it. The
/// current it deposited moved only part of its charge off the boundary nodes;
/// the rest would stay there as virtual charge, breaking Gauss's law. This
/// continues the particle's path to the ghost cell's outer face, normal to each
/// open boundary it is beyond (one after the other, for corners), so all of its
/// charge leaves the domain's nodes.
///
/// Injects nothing, but must run where injectors do: after the particle push
/// and before particle boundary exchange.
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`
template <typename PUSH_PARTICLES>
class OpenBoundaryOutflow
  : public InjectorBase<typename PUSH_PARTICLES::Mparticles,
                        typename PUSH_PARTICLES::MfieldsState>
{
public:
  using PushParticles = PUSH_PARTICLES;
  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using Current = typename PushParticles::Current;
  using real_t = typename PushParticles::real_t;
  using Real3 = Vec3<real_t>;

  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    const Grid_t& grid = mprts.grid();
    Real3 dxi{real_t(grid.domain.dx_inv[0]), real_t(grid.domain.dx_inv[1]),
              real_t(grid.domain.dx_inv[2])};
    Current current(grid);
    auto accessor = mprts.accessor();

    for (int p = 0; p < grid.n_patches(); p++) {
      bool open_lo[3], open_hi[3];
      bool any_open = false;
      for (int d = 0; d < 3; d++) {
        open_lo[d] =
          grid.atBoundaryLo(p, d) && grid.bc.prt_lo[d] == BND_PRT_OPEN;
        open_hi[d] =
          grid.atBoundaryHi(p, d) && grid.bc.prt_hi[d] == BND_PRT_OPEN;
        any_open |= open_lo[d] || open_hi[d];
      }
      if (!any_open) {
        continue;
      }

      auto flds = mflds[p];
      typename Current::fields_t J(flds);

      for (auto prt : accessor[p]) {
        Real3 pos = Real3(prt.x()) * dxi;
        for (int d = 0; d < 3; d++) {
          real_t outer_face;
          if (open_lo[d] && pos[d] < 0) {
            outer_face = -1;
          } else if (open_hi[d] && pos[d] >= grid.ldims[d]) {
            outer_face = grid.ldims[d] + 1;
          } else {
            continue;
          }

          // the cell the particle is in; an earlier move may have left it on
          // the outer face of an upper ghost cell, which belongs to that cell
          Int3 idx = pos.fint();
          for (int e = 0; e < 3; e++) {
            if (open_hi[e]) {
              idx[e] = std::min(idx[e], grid.ldims[e]);
            }
          }

          Real3 initial_pos = pos;
          Real3 final_pos = pos.with_component(d, outer_face);
          Real3 no_v{0.f, 0.f, 0.f};
          current.calc_j(J, initial_pos, final_pos, idx, idx,
                         real_t(prt.q() * prt.w()), no_v);
          pos[d] = outer_face;
        }
      }
    }
  }
};
