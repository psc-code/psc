#pragma once

#include "psc.h"
#include "../axis.hxx"
#include "kg/Vec3.h"
#include "field_bc_base.hxx"
#include "field_bc_util.hxx"
#include "../psc_bnd/psc_bnd_util.hxx"

namespace psc
{
namespace bnd
{
namespace field
{

template <typename Dim, typename MfieldsState>
struct ConductingWall : FieldBcBase<MfieldsState>
{
  using dim_t = Dim;

  ConductingWall(Axis d, LoHi lohi) : d{d}, lohi{lohi} {}

  void apply_j_bcs(MfieldsState& mflds) override
  {
    // todo
  }

  void apply_e_bcs(MfieldsState& mflds) override
  {
    const auto& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      if (lohi == Lo && grid.atBoundaryLo(p, d)) {
        detail::set_lower_ghosts_to_nan<dim_t>(mflds, p, d, EX, true);

        auto F = make_Fields3d<Dim>(mflds[p]);
        const int* ldims = mflds.grid().ldims;
        Int3 ib = mflds.ib(), im = mflds.im();

        if (d == Axis::Y) {
          for (int iz = -2; iz < ldims[2] + 2; iz++) {
            // FIXME, needs to be for other dir, too, and it's ugly
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(EX, ix, 0, iz) = 0.;
              F(EX, ix, -1, iz) = F(EX, ix, 1, iz);

              F(EY, ix, -1, iz) = -F(EY, ix, 0, iz);

              F(EZ, ix, 0, iz) = 0.;
              F(EZ, ix, -1, iz) = F(EZ, ix, 1, iz);
            }
          }
        } else if (d == Axis::Z) {
          for (int iy = -2; iy < ldims[1] + 2; iy++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(EX, ix, iy, 0) = 0.;
              F(EX, ix, iy, -1) = F(EX, ix, iy, 1);

              F(EY, ix, iy, 0) = 0.;
              F(EY, ix, iy, -1) = F(EY, ix, iy, 1);

              F(EZ, ix, iy, -1) = -F(EZ, ix, iy, 0);
            }
          }
        } else {
          assert(0);
        }
      }

      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        detail::set_upper_ghosts_to_nan<dim_t>(mflds, p, d, EX, true);

        auto F = make_Fields3d<Dim>(mflds[p]);
        const int* ldims = mflds.grid().ldims;
        Int3 ib = mflds.ib(), im = mflds.im();

        if (d == Axis::Y) {
          int my _mrc_unused = ldims[1];
          for (int iz = -2; iz < ldims[2] + 2; iz++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(EX, ix, my, iz) = 0.;
              F(EX, ix, my + 1, iz) = F(EX, ix, my - 1, iz);

              F(EY, ix, my, iz) = -F(EY, ix, my - 1, iz);

              F(EZ, ix, my, iz) = 0.;
              F(EZ, ix, my + 1, iz) = F(EZ, ix, my - 1, iz);
            }
          }
        } else if (d == Axis::Z) {
          int mz = ldims[2];
          for (int iy = -2; iy < ldims[1] + 2; iy++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(EX, ix, iy, mz) = 0.;
              F(EX, ix, iy, mz + 1) = F(EX, ix, iy, mz - 1);

              F(EY, ix, iy, mz) = 0.;
              F(EY, ix, iy, mz + 1) = F(EY, ix, iy, mz - 1);

              F(EZ, ix, iy, mz) = -F(EZ, ix, iy, mz - 1);
            }
          }
        } else {
          assert(0);
        }
      }
    }
  }

  void apply_h_bcs(MfieldsState& mflds) override
  {
    // todo
  }

  Axis d;
  LoHi lohi;
};

} // namespace field
} // namespace bnd
} // namespace psc
