#pragma once

#include "psc.h"
#include "../axis.hxx"
#include "kg/Vec3.h"
#include "kg/VecRange.hxx"
#include "field_bc_base.hxx"
#include "field_bc_util.hxx"
#include "../psc_bnd/psc_bnd_util.hxx"

namespace psc
{
namespace bnd
{
namespace field
{

template <typename Dim, typename MfieldsState, typename P>
struct Radiating : FieldBcBase<MfieldsState>
{
  using dim_t = Dim;
  using Pulse = P;
  using real_t = typename MfieldsState::real_t;
  using Real3 = typename MfieldsState::Real3;

  Radiating(Pulse pulse, Axis d, LoHi lohi) : pulse{pulse}, d{d}, lohi{lohi} {}

  void apply_j_bcs(MfieldsState& mflds) override {}

  void apply_e_bcs(MfieldsState& mflds) override {}

  void apply_h_bcs(MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();

    pulse.tick(grid.time());

    for (int p = 0; p < mflds.n_patches(); p++) {
      if (lohi == Lo && grid.atBoundaryLo(p, d)) {
        psc::bnd::field::detail::set_lower_ghosts_to_nan<dim_t>(mflds, p, d, HX,
                                                                false);

        auto F = make_Fields3d<dim_t>(mflds[p]);
        const Grid_t& grid = mflds.grid();
        real_t dt = grid.dt;
        Real3 dtdx = dt * Real3(grid.domain.dx_inv);

        int d0 = (d + 0) % 3, d1 = (d + 1) % 3, d2 = (d + 2) % 3;
        int H0 = HX + d0, H1 = HX + d1, H2 = HX + d2;
        int E1 = EX + d1, E2 = EX + d2;
        int J1 = JXI + d1, J2 = JXI + d2;

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d0] = -1;
        stop[d0] = 0;

        for (Int3 i3 : VecRange(start, stop)) {
          Int3 edge_idx = i3 + Int3::unit(d0);

          real_t s = 0.0;
          real_t p = 0.0;

          Real3 x3_s = (Real3(edge_idx) + Real3::unit(d1) * real_t(0.5)) *
                       Real3(grid.domain.dx);
          Real3 x3_p = (Real3(edge_idx) + Real3::unit(d2) * real_t(0.5)) *
                       Real3(grid.domain.dx);

          s = pulse.sample_exterior_field_lo(EX + d1, grid.time(), p, x3_s) +
              pulse.sample_exterior_field_lo(HX + d2, grid.time(), p, x3_s);
          p = pulse.sample_exterior_field_lo(EX + d2, grid.time(), p, x3_p) -
              pulse.sample_exterior_field_lo(HX + d1, grid.time(), p, x3_p);

          F(H2, i3) =
            (2.f * s - 2.f * F(E1, edge_idx) -
             dtdx[d2] * (F(H0, edge_idx) - F(H0, edge_idx - Int3::unit(d2))) -
             (1.f - dtdx[d0]) * F(H2, edge_idx) + dt * F(J1, edge_idx)) /
            (1.f + dtdx[d0]);
          F(H1, i3) =
            (-2.f * p + 2.f * F(E2, edge_idx) -
             dtdx[d1] * (F(H0, edge_idx) - F(H0, edge_idx - Int3::unit(d1))) -
             (1.f - dtdx[d0]) * F(H1, edge_idx) - dt * F(J2, edge_idx)) /
            (1.f + dtdx[d0]);
        }
      }

      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        psc::bnd::field::detail::set_upper_ghosts_to_nan<dim_t>(mflds, p, d, HX,
                                                                false);

        auto F = make_Fields3d<dim_t>(mflds[p]);
        const Grid_t& grid = mflds.grid();
        Int3 ldims = grid.ldims;
        real_t dt = grid.dt;
        Real3 dtdx = dt * Real3(grid.domain.dx_inv);

        int d0 = (d + 0) % 3, d1 = (d + 1) % 3, d2 = (d + 2) % 3;
        int H0 = HX + d0, H1 = HX + d1, H2 = HX + d2;
        int E1 = EX + d1, E2 = EX + d2;
        int J1 = JXI + d1, J2 = JXI + d2;

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d0] = grid.ldims[d0];
        stop[d0] = grid.ldims[d0] + 1;

        for (Int3 i3 : VecRange(start, stop)) {
          Int3 edge_idx = i3 - Int3::unit(d0);

          real_t s = 0.0;
          real_t p = 0.0;

          Real3 x3_s = (Real3(edge_idx) + Real3::unit(d1) * real_t(0.5)) *
                       Real3(grid.domain.dx);
          Real3 x3_p = (Real3(edge_idx) + Real3::unit(d2) * real_t(0.5)) *
                       Real3(grid.domain.dx);

          s = pulse.sample_exterior_field_hi(EX + d1, grid.time(), p, x3_s) -
              pulse.sample_exterior_field_hi(HX + d2, grid.time(), p, x3_s);
          p = pulse.sample_exterior_field_hi(EX + d2, grid.time(), p, x3_p) +
              pulse.sample_exterior_field_hi(HX + d1, grid.time(), p, x3_p);

          F(H2, i3) = (-2.f * s + 2.f * F(E1, i3) +
                       dtdx[d2] * (F(H0, i3) - F(H0, i3 - Int3::unit(d2))) -
                       (1.f - dtdx[d0]) * F(H2, edge_idx) - dt * F(J1, i3)) /
                      (1.f + dtdx[d0]);
          F(H1, i3) = (2.f * p - 2.f * F(E2, i3) +
                       dtdx[d1] * (F(H0, i3) - F(H0, i3 - Int3::unit(d1))) -
                       (1.f - dtdx[d0]) * F(H1, edge_idx) + dt * F(J2, i3)) /
                      (1.f + dtdx[d0]);
        }
      }
    }
  }

  Axis d;
  LoHi lohi;

private:
  Pulse pulse;
};

} // namespace field
} // namespace bnd
} // namespace psc
