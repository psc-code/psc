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

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d] = 0;
        stop[d] = start[d] + 1;

        int E0 = EX + d;
        int E1 = EX + d.next();
        int E2 = EX + d.prev();

        for (Int3 i3 : VecRange(start, stop)) {
          // 1. transverse components: wall is at relative index 0
          F(E1, i3) = 0.0;
          F(E2, i3) = 0.0;
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(E1, i3 - j3) = F(E1, i3 + j3);
            F(E2, i3 - j3) = F(E2, i3 + j3);
          }

          // 2. normal component: wall is at relative index -1/2
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(E0, i3 - j3) = -F(E0, i3 + j3 - Int3::unit(d));
          }
        }
      }

      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        detail::set_upper_ghosts_to_nan<dim_t>(mflds, p, d, EX, true);

        auto F = make_Fields3d<Dim>(mflds[p]);
        Int3 ldims = mflds.grid().ldims;

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d] = ldims[d];
        stop[d] = start[d] + 1;

        int E0 = EX + d;
        int E1 = EX + d.next();
        int E2 = EX + d.prev();

        for (Int3 i3 : VecRange(start, stop)) {
          // 1. transverse components: wall is at relative index n=ldims[d]
          F(E1, i3) = 0.0;
          F(E2, i3) = 0.0;
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(E1, i3 + j3) = F(E1, i3 - j3);
            F(E2, i3 + j3) = F(E2, i3 - j3);
          }

          // 2. normal component: wall is at relative index n-1/2
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(E0, i3 + j3) = -F(E0, i3 - j3 + Int3::unit(d));
          }
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
