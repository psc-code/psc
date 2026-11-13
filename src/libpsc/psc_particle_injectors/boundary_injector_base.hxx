#pragma once

#include "grid.hxx"
#include "particle.h"
#include "injector_base.hxx"
#include "../axis.hxx"
#include "../bits.hxx"
#include "../psc_bnd/psc_bnd_util.hxx"

namespace psc
{
namespace bnd
{
namespace particle
{

/// normalized coordinate along `d` of the first ghost cell's face away from the
/// domain, at the given boundary
inline double outer_ghost_face(const Grid_t& grid, int d, LoHi lohi)
{
  return lohi == LoHi::Lo ? -1 : grid.ldims[d] + 1;
}

/// @brief Common base of injectors that inject particles through one boundary,
/// normal to axis `d`, at its lower or upper end.
/// @tparam PUSH_PARTICLES type that provides the types `Mparticles`,
/// `MfieldsState`, `Current`, `real_t`
template <typename PUSH_PARTICLES>
class BoundaryInjectorBase
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

  BoundaryInjectorBase(Axis d, LoHi lohi) : d{d}, lohi{lohi} {}

protected:
  bool lo() const { return lohi == LoHi::Lo; }

  /// +1 or -1, whichever points from the boundary into the domain
  int inwards() const { return lo() ? 1 : -1; }

  bool at_boundary(const Grid_t& grid, int p) const
  {
    return lo() ? grid.atBoundaryLo(p, d) : grid.atBoundaryHi(p, d);
  }

  /// index along `d` of the cells adjacent to the boundary
  int edge_idx(const Grid_t& grid) const
  {
    return lo() ? 0 : grid.ldims[d] - 1;
  }

  void check_axis(const Grid_t& grid) const
  {
    if (grid.isInvar(d)) {
      LOG_ERROR("can't inject through boundaries of invariant dimension %d\n",
                int(d));
    }
  }

  /// Injects `prt`, which moved with velocity `v` from
  /// `initial_normalized_pos`, in ghost cell `initial_idx`, to its current
  /// position, and deposits its current as if it came from the ghost cell's
  /// outer face. Its virtual charge then lies entirely beyond the boundary
  /// node, so Gauss's law holds there.
  template <typename Injector>
  void inject_with_current(const Grid_t& grid, Injector& injector,
                           Current& current, typename Current::fields_t& J,
                           const psc::particle::Inject& prt,
                           Real3 initial_normalized_pos, Int3 initial_idx,
                           Real3 v, real_t charge) const
  {
    injector.inject_local_unchecked(prt);

    Real3 dxi = Real3(grid.domain.dx_inv);
    Real3 final_normalized_pos = Real3(prt.x) * dxi;
    Int3 final_idx = final_normalized_pos.fint();
    initial_normalized_pos[d] = outer_ghost_face(grid, d, lohi);

    current.calc_j(J, initial_normalized_pos, final_normalized_pos, final_idx,
                   initial_idx, charge, v);
  }

  Axis d;
  LoHi lohi;
};

} // namespace particle
} // namespace bnd
} // namespace psc
