
#pragma once

// T. Zh. Esirkepov, "Exact charge conservation scheme for particle-in-cell
// simulation with an arbitrary form-factor", Computer Physics Communications
// 135 (2001) 144

#include "pushp_current_esirkepov.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"

namespace
{
/**
 * @brief Handle exiting particles. When a particle crosses an open boundary, it
 * is dropped. Before that happens, an exiting current needs to be deposited, or
 * else virtual charge will accumulate at the boundary.
 *
 * This function currently asserts that there are no open boundaries (and is
 * thus a noop in release builds). See `push_particles_1vb.hxx` for an actual
 * implementation for 1st-order particles.
 * @tparam Real 'float' or 'double'
 * @param grid the grid
 * @param p the patch index
 */
template <typename Real>
void handle_exiting(const Grid_<Real>& grid, int p)
{
  for (int d = 0; d < 3; d++) {
    assert(grid.bc.prt_lo[d] != BND_PRT_OPEN);
    assert(grid.bc.prt_hi[d] != BND_PRT_OPEN);
  }
}
} // namespace

// ======================================================================
// PushParticlesEsirkepov

template <typename C>
struct PushParticlesEsirkepov
{
  static const int MAX_NR_KINDS = 10;

  using Mparticles = typename C::Mparticles;
  using MfieldsState = typename C::MfieldsState;
  using AdvanceParticle_t = typename C::AdvanceParticle_t;
  using InterpolateEM_t = typename C::InterpolateEM_t;
  using Current =
    CurrentEsirkepov<typename C::Order, typename C::Dim,
                     Fields3d<typename MfieldsState::fields_view_t::Storage>,
                     InterpolateEM_t>;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;

  using checks_order =
    checks_order_2nd; // FIXME, sometimes 1st even with Esirkepov

  // ----------------------------------------------------------------------
  // push_mprts

  static void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    Real3 dxi = grid.domain.dx_inv;
    real_t dq_kind[MAX_NR_KINDS];
    auto& kinds = grid.kinds;
    assert(kinds.size() <= MAX_NR_KINDS);
    for (int k = 0; k < kinds.size(); k++) {
      dq_kind[k] = .5f * grid.norm.eta * grid.dt * kinds[k].q / kinds[k].m;
    }
    InterpolateEM_t ip;
    AdvanceParticle_t advance(grid.dt);
    Current current(grid);

    auto accessor = mprts.accessor_();
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = mflds[p];
      auto prts = accessor[p];
      typename InterpolateEM_t::fields_t EM(flds.storage(), flds.ib());
      typename Current::fields_t J(flds.storage(), flds.ib());

      flds.storage().view(_all, _all, _all, _s(JXI, JXI + 3)) = real_t(0);

      for (auto prt : prts) {
        Real3& x = prt.x();

        real_t xm[3];
        for (int d = 0; d < 3; d++) {
          xm[d] = x[d] * dxi[d];
        }
        ip.set_coeffs(xm);

        // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt
        current.charge_before(ip);

        // FIELD INTERPOLATION
        Real3 E = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
        Real3 H = {ip.hx(EM), ip.hy(EM), ip.hz(EM)};

        // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
        real_t dq = dq_kind[prt.kind()];
        advance.push_p(prt.u(), E, H, dq);

        // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
        auto v = advance.calc_v(prt.u());
        advance.push_x(x, v);

        // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt
        current.charge_after(x);

        // CURRENT DENSITY AT (n+1.0)*dt
        current.prep(prt.qni_wni(), v);
        current.calc(J);

        handle_exiting(grid, p);
      }
    }
  }
};
