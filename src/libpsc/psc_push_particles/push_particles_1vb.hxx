
// J. Villasenor and O. Buneman, "Rigorous charge conservation for local
// electromagnetic field solvers", Computer Physics Communications 69 (1992) 306

namespace
{
/**
 * @brief "Pushes" an exiting particle to the outer edge of the first layer of
 * ghost cells. This should happen before current deposition, so the full
 * exiting current is deposited. Note that higher order particles would need to
 * be pushed deeper into the ghost region.
 *
 * Actually, this function works with the grid-normalized particle position, not
 * the particle itself (or even its actual position). The grid-normalized
 * position is what's actually used to deposit current, and the particle is
 * dropped later, so its position doesn't need to be updated.
 * @tparam Real `float` or `double`
 * @param final_x_normed the particle's grid-normalized position, which is
 * possibly mutated
 * @param final_i3 the particle's 3d cell index
 * @param grid the grid
 * @param p the patch index
 */
template <typename Real>
void exit_to_edge(Vec3<Real>& final_x_normed, const Int3& final_i3,
                  const Grid_t& grid, int p)
{
  // FIXME #948112531345 (also see the other FIXMEs with this id)
  // Current deposited in ghost corners isn't sent to other patches, and
  // thus leads to false-positive gauss errors.
  // This is a general problem with non-periodic boundaries that have currents
  // in ghost cells (i.e., just open BCs as of now), not just this function.
  for (int d = 0; d < 3; d++) {
    if (grid.bc.prt_lo[d] == BND_PRT_OPEN && grid.atBoundaryLo(p, d)) {
      if (final_i3[d] < 0) {
        final_x_normed[d] = -1.0;
      }
    }

    if (grid.bc.prt_hi[d] == BND_PRT_OPEN && grid.atBoundaryHi(p, d)) {
      if (final_i3[d] >= grid.ldims[d]) {
        final_x_normed[d] = grid.ldims[d] + 1.0;
      }
    }
  }
}
} // namespace

// ======================================================================
// PushParticlesVb

template <typename C>
struct PushParticlesVb
{
  static const int MAX_NR_KINDS = 10;

  using Mparticles = typename C::Mparticles;
  using MfieldsState = typename C::MfieldsState;
  using AdvanceParticle_t = typename C::AdvanceParticle_t;
  using InterpolateEM_t = typename C::InterpolateEM_t;
  using Current = typename C::Current_t;
  using Dim = typename C::Dim;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;

  using checks_order = checks_order_1st;

  // ----------------------------------------------------------------------
  // push_mprts

  static void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    Real3 dxi = Real3(grid.domain.dx).inv();
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
      typename Current::fields_t J(flds);

      flds.storage().view(_all, _all, _all, _s(JXI, JXI + 3)) = real_t(0);

      for (auto prt : prts) {
        Real3 initial_pos_normalized = prt.x() * dxi;
        ip.set_coeffs(initial_pos_normalized);

        // FIELD INTERPOLATION
        Real3 E = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
        Real3 H = {ip.hx(EM), ip.hy(EM), ip.hz(EM)};

        // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
        real_t dq = dq_kind[prt.kind()];
        advance.push_p(prt.u(), E, H, dq);

        // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
        auto v = advance.calc_v(prt.u());
        advance.push_x(prt.x(), v);

        Real3 final_pos_normalized = prt.x() * dxi;
        Int3 final_index = final_pos_normalized.fint();

        exit_to_edge(final_pos_normalized, final_index, grid, p);

        // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
        Int3 initial_index;
        if (!Dim::InvarX::value) {
          initial_index[0] = ip.cx.g.l;
        }
        if (!Dim::InvarY::value) {
          initial_index[1] = ip.cy.g.l;
        }
        if (!Dim::InvarZ::value) {
          initial_index[2] = ip.cz.g.l;
        }
        current.calc_j(J, initial_pos_normalized, final_pos_normalized,
                       final_index, initial_index, prt.qni_wni(), v);
      }
    }
  }

  // ----------------------------------------------------------------------
  // stagger_mprts_patch

  static void stagger_mprts_patch(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    Real3 dxi = grid.domain.dx_inv;
    real_t dq_kind[MAX_NR_KINDS];
    auto& kinds = grid.kinds;
    assert(kinds.size() <= MAX_NR_KINDS);
    for (int k = 0; k < kinds.size(); k++) {
      dq_kind[k] = .5f * grid.eta * grid.dt * kinds[k].q / kinds[k].m;
    }
    InterpolateEM_t ip;
    AdvanceParticle_t advance(grid.dt);

    auto accessor = mprts.accessor_();
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = mflds[p];
      auto prts = accessor[p];
      typename InterpolateEM_t::fields_t EM(flds);

      for (auto prt : prts) {
        // field interpolation
        real_t* x = prt.x;

        real_t initial_pos_normalized[3];
        for (int d = 0; d < 3; d++) {
          initial_pos_normalized[d] = x[d] * dxi[d];
        }

        // FIELD INTERPOLATION

        ip.set_coeffs(initial_pos_normalized);
        // FIXME, we're not using EM instead flds_em
        real_t E[3] = {ip.ex(EM), ip.ey(EM), ip.ez(EM)};
        real_t H[3] = {ip.hx(EM), ip.hy(EM), ip.hz(EM)};

        // x^(n+1/2), p^{n+1/2} -> x^(n+1/2), p^{n}
        int kind = prt->kind();
        real_t dq = dq_kind[kind];
        advance.push_p(&prt->pxi, E, H, -.5f * dq);
      }
    }
  }
};
