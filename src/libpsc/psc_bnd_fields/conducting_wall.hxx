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

/**
 * @brief A perfect electrical conductor. This implementation assumes particles
 * are specularly reflected.
 * @tparam Dim dimension type
 * @tparam MfieldsState fields type
 */
template <typename Dim, typename MfieldsState>
struct ConductingWall : FieldBcBase<MfieldsState>
{
  using dim_t = Dim;

  ConductingWall(Axis d, LoHi lohi) : d{d}, lohi{lohi} {}

  void apply_j_bcs(MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      if (lohi == Lo && grid.atBoundaryLo(p, d)) {
        auto F = make_Fields3d<dim_t>(mflds[p]);
        const int* ldims = mflds.grid().ldims;
        Int3 ib = mflds.ib(), im = mflds.im();

        if (d == 1) {
          for (int iz = -2; iz < ldims[2] + 2; iz++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(JXI, ix, 1, iz) += F(JXI, ix, -1, iz);
              F(JXI, ix, -1, iz) = 0.;

              F(JYI, ix, 0, iz) -= F(JYI, ix, -1, iz);
              F(JYI, ix, -1, iz) = 0.;

              F(JZI, ix, 1, iz) += F(JZI, ix, -1, iz);
              F(JZI, ix, -1, iz) = 0.;
            }
          }
        } else if (d == 2) {
          for (int iy = -2; iy < ldims[1] + 2; iy++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(JXI, ix, iy, 1) += F(JXI, ix, iy, -1);
              F(JXI, ix, iy, -1) = 0.;

              F(JYI, ix, iy, 1) += F(JYI, ix, iy, -1);
              F(JYI, ix, iy, -1) = 0.;

              F(JZI, ix, iy, 0) -= F(JZI, ix, iy, -1);
              F(JZI, ix, iy, -1) = 0.;
            }
          }
        } else {
          assert(0);
        }
      }

      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        auto F = make_Fields3d<dim_t>(mflds[p]);
        const int* ldims = mflds.grid().ldims;
        Int3 ib = mflds.ib(), im = mflds.im();

        if (d == 1) {
          int my _mrc_unused = ldims[1];
          for (int iz = -2; iz < ldims[2] + 2; iz++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(JXI, ix, my - 1, iz) += F(JXI, ix, my + 1, iz);
              F(JXI, ix, my + 1, iz) = 0.;

              F(JYI, ix, my - 1, iz) -= F(JYI, ix, my, iz);
              F(JYI, ix, my, iz) = 0.;

              F(JZI, ix, my - 1, iz) += F(JZI, ix, my + 1, iz);
              F(JZI, ix, my + 1, iz) = 0.;
            }
          }
        } else if (d == 2) {
          int mz = ldims[2];
          for (int iy = -2; iy < ldims[1] + 2; iy++) {
            for (int ix = std::max(-2, ib[0]);
                 ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
              F(JXI, ix, iy, mz - 1) += F(JXI, ix, iy, mz + 1);
              F(JXI, ix, iy, mz + 1) = 0.;

              F(JYI, ix, iy, mz - 1) += F(JYI, ix, iy, mz + 1);
              F(JYI, ix, iy, mz + 1) = 0.;

              F(JZI, ix, iy, mz - 1) -= F(JZI, ix, iy, mz);
              F(JZI, ix, iy, mz) = 0.;
            }
          }
        } else {
          assert(0);
        }
      }
    }
  }

  /**
   * @brief Set transverse E at the wall's surface to 0. Nominally, E would be 0
   * at interior points too, but instead we mirror E so that reflecting
   * particles feel the "right" electric force.
   *
   * Note that the normal force is reflected in the interior. To see why this is
   * necessary, consider a particle with v=0 at the surface. If the normal E
   * force on it is nonzero and towards the wall, the particle will
   * spontaneously bounce away from the wall. It will then return to the
   * wall—with more kinetic energy than before—and bounce again. Flipping normal
   * E avoids this runaway effect.
   * @param mflds fields
   */
  void apply_e_bcs(MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();

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

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d] = grid.ldims[d];
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

  /**
   * @brief Normal H at the wall's surface is a no-op, but set interior H such
   * that reflecting particles feel the "right" magnetic force.
   *
   * Transverse H is flipped to ensure that a particle at the wall's surface
   * with nominally-nonzero normal velocity—and thus, actually zero average
   * normal velocity, since half of the cloud is reflected—experiences no
   * magnetic force.
   * @param mflds fields
   */
  void apply_h_bcs(MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      if (lohi == Lo && grid.atBoundaryLo(p, d)) {
        detail::set_lower_ghosts_to_nan<dim_t>(mflds, p, d, HX, false);

        auto F = make_Fields3d<dim_t>(mflds[p]);

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d] = 0;
        stop[d] = start[d] + 1;

        int H0 = HX + d;
        int H1 = HX + d.next();
        int H2 = HX + d.prev();

        for (Int3 i3 : VecRange(start, stop)) {
          // 1. transverse components: wall is at relative index -1/2
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(H1, i3 - j3) = -F(H1, i3 + j3 - Int3::unit(d));
            F(H2, i3 - j3) = -F(H2, i3 + j3 - Int3::unit(d));
          }

          // 2. normal component: wall is at relative index 0
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(H0, i3 - j3) = F(H0, i3 + j3);
          }
        }
      }

      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        detail::set_upper_ghosts_to_nan<dim_t>(mflds, p, d, HX, false);

        auto F = make_Fields3d<dim_t>(mflds[p]);

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d] = grid.ldims[d];
        stop[d] = start[d] + 1;

        int H0 = HX + d;
        int H1 = HX + d.next();
        int H2 = HX + d.prev();

        for (Int3 i3 : VecRange(start, stop)) {
          // 1. transverse components: wall is at relative index n-1/2
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(H1, i3 + j3) = -F(H1, i3 - j3 + Int3::unit(d));
            F(H2, i3 + j3) = -F(H2, i3 - j3 + Int3::unit(d));
          }

          // 2. normal component: wall is at relative index n=ldims[d]
          for (Int3 j3 = Int3::unit(d); j3[d] <= mflds.ibn()[d]; j3[d]++) {
            F(H0, i3 + j3) = F(H0, i3 - j3);
          }
        }
      }
    }
  }

  Axis d;
  LoHi lohi;
};

} // namespace field
} // namespace bnd
} // namespace psc
