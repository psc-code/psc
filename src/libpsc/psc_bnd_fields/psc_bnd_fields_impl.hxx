
#include "psc.h"
#include "kg/VecRange.hxx"
#include "fields.hxx"
#include "bnd_fields.hxx"
#include "radiating_bnd.hxx"

#include <mrc_bits.h>

#include <limits>

// #define DEBUG

template <typename MFIELDS_STATE, typename Dim>
struct BndFields_ : BndFieldsBase
{
  using Self = BndFields_<MFIELDS_STATE, Dim>;
  using MfieldsState = MFIELDS_STATE;
  using real_t = typename MfieldsState::real_t;
  using Real3 = Vec3<real_t>;
  using fields_view_t = typename MfieldsState::fields_view_t;
  using dim_t = Dim;

  // ----------------------------------------------------------------------
  // fill_ghosts_E

  void fill_ghosts_E(MfieldsState& mflds)
  {
    const auto& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      // lo
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryLo(p, d)) {
          switch (grid.bc.fld_lo[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              conducting_wall_E_lo(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              set_lower_ghosts(mflds, p, d, EX, background_e, false);
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }

      // hi
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryHi(p, d)) {
          switch (grid.bc.fld_hi[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              conducting_wall_E_hi(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              set_upper_ghosts(mflds, p, d, EX, background_e, false);
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // fill_ghosts_H

  void fill_ghosts_H(MfieldsState& mflds)
  {
    const auto& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      // lo
      for (int d = 0; d < 3; d++) {
        if (grid.bc.fld_lo[d] == BND_FLD_OPEN && radiation) {
          radiation->update_cache_lower(grid.time(), d);
        }

        if (grid.atBoundaryLo(p, d)) {
          switch (grid.bc.fld_lo[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              conducting_wall_H_lo(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              radiative_H_lo(mflds, p, d);
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
      // hi
      for (int d = 0; d < 3; d++) {
        if (grid.bc.fld_hi[d] == BND_FLD_OPEN && radiation) {
          radiation->update_cache_upper(grid.time(), d);
        }

        if (grid.atBoundaryHi(p, d)) {
          switch (grid.bc.fld_hi[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              conducting_wall_H_hi(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              radiative_H_hi(mflds, p, d);
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // add_ghosts_J

  void add_ghosts_J(MfieldsState& mflds)
  {
    const auto& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      // lo
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryLo(p, d)) {
          switch (grid.bc.fld_lo[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              conducting_wall_J_lo(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
      // hi
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryHi(p, d)) {
          switch (grid.bc.fld_hi[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              conducting_wall_J_hi(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
    }
  }

  static void set_lower_ghosts_to_nan(MfieldsState& mflds, int p, int d, int mb,
                                      bool include_edge)
  {
#ifndef DEBUG
    return;
#endif
    real_t nan = std::numeric_limits<real_t>::quiet_NaN();
    set_lower_ghosts(mflds, p, d, mb, {nan, nan, nan}, include_edge);
  }

  static void set_upper_ghosts_to_nan(MfieldsState& mflds, int p, int d, int mb,
                                      bool include_edge)
  {
#ifndef DEBUG
    return;
#endif
    real_t nan = std::numeric_limits<real_t>::quiet_NaN();
    set_upper_ghosts(mflds, p, d, mb, {nan, nan, nan}, include_edge);
  }

  /**
   * @brief Set E or B lower ghosts to the given constants (each component has
   * its own constant).
   * @param mflds mflds
   * @param p patch index
   * @param d which dimension to set the ghosts of
   * @param mb `EX` or `HX`; note that `mb+1` and `mb+2` are also set
   * @param val the constants
   * @param include_edge whether or not values located on exact domain edges
   * should be considered "ghosts"
   */
  static void set_lower_ghosts(MfieldsState& mflds, int p, int d, int mb,
                               Real3 val, bool include_edge)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    Int3 start = mflds.ib();
    Int3 stop = mflds.ib() + mflds.im();
    stop[d] = 0;

    // TODO use gtensor views instead of VecRange

    for (int m = mb; m < mb + 3; m++) {
      for (Int3 i3 : VecRange(start, stop)) {
        F(m, i3) = val[m - mb];
      }
    }

    if (!include_edge) {
      return;
    }

    Int3 edge_start = mflds.ib();
    Int3 edge_stop = mflds.ib() + mflds.im();
    edge_start[d] = 0;
    edge_stop[d] = 1;

    for (int m = mb; m < mb + 3; m++) {
      bool edge_ec = mb == EX && m - mb != d;
      bool edge_fc = mb == HX && m - mb == d;

      if (edge_ec || edge_fc) {
        for (Int3 i3 : VecRange(edge_start, edge_stop)) {
          F(m, i3) = val[m - mb];
        }
      }
    }
  }

  /**
   * @brief Set E or B upper ghosts to the given constants (each component has
   * its own constant).
   * @param mflds mflds
   * @param p patch index
   * @param d which dimension to set the ghosts of
   * @param mb `EX` or `HX`; note that `mb+1` and `mb+2` are also set
   * @param val the constants
   * @param include_edge whether or not values located on exact domain edges
   * should be considered "ghosts"
   */
  static void set_upper_ghosts(MfieldsState& mflds, int p, int d, int mb,
                               Real3 val, bool include_edge)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    Int3 start = mflds.ib();
    Int3 stop = mflds.ib() + mflds.im();
    start[d] = mflds.grid().ldims[d] + 1;

    // TODO use gtensor views instead of VecRange

    for (int m = mb; m < mb + 3; m++) {
      for (Int3 i3 : VecRange(start, stop)) {
        F(m, i3) = val[m - mb];
      }
    }

    Int3 edge_start = mflds.ib();
    Int3 edge_stop = mflds.ib() + mflds.im();
    edge_start[d] = mflds.grid().ldims[d];
    edge_stop[d] = mflds.grid().ldims[d] + 1;

    for (int m = mb; m < mb + 3; m++) {
      bool not_edge_ec = mb == EX && m - mb == d;
      bool not_edge_fc = mb == HX && m - mb != d;

      if (not_edge_ec || not_edge_fc || include_edge)
        for (Int3 i3 : VecRange(edge_start, edge_stop)) {
          F(m, i3) = val[m - mb];
        }
    }
  }

  void conducting_wall_E_lo(MfieldsState& mflds, int p, int d)
  {
    set_lower_ghosts_to_nan(mflds, p, d, EX, true);

    auto F = make_Fields3d<dim_t>(mflds[p]);
    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
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
    } else if (d == 2) {
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

  void conducting_wall_E_hi(MfieldsState& mflds, int p, int d)
  {
    set_upper_ghosts_to_nan(mflds, p, d, EX, true);

    auto F = make_Fields3d<dim_t>(mflds[p]);
    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
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
    } else if (d == 2) {
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

  void conducting_wall_H_lo(MfieldsState& mflds, int p, int d)
  {
    set_lower_ghosts_to_nan(mflds, p, d, HX, false);

    auto F = make_Fields3d<dim_t>(mflds[p]);
    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, -1, iz) = -F(HX, ix, 0, iz);

          F(HY, ix, -1, iz) = F(HY, ix, 1, iz);

          F(HZ, ix, -1, iz) = -F(HZ, ix, 0, iz);
        }
      }
    } else if (d == 2) {
      for (int iy = -2; iy < ldims[1] + 2; iy++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, iy, -1) = -F(HX, ix, iy, 0);

          F(HY, ix, iy, -1) = -F(HY, ix, iy, 0);

          F(HZ, ix, iy, -1) = F(HZ, ix, iy, 1);
        }
      }
    } else {
      assert(0);
    }
  }

  void conducting_wall_H_hi(MfieldsState& mflds, int p, int d)
  {
    set_upper_ghosts_to_nan(mflds, p, d, HX, false);

    auto F = make_Fields3d<dim_t>(mflds[p]);

    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
      int my _mrc_unused = ldims[1];
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, my, iz) = -F(HX, ix, my - 1, iz);

          F(HY, ix, my + 1, iz) = F(HY, ix, my - 1, iz);

          F(HZ, ix, my, iz) = -F(HZ, ix, my - 1, iz);
        }
      }
    } else if (d == 2) {
      int mz = ldims[2];
      for (int iy = -2; iy < ldims[1] + 2; iy++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, iy, mz) = -F(HX, ix, iy, mz - 1);

          F(HY, ix, iy, mz) = -F(HY, ix, iy, mz - 1);

          F(HZ, ix, iy, mz + 1) = F(HZ, ix, iy, mz - 1);
        }
      }
    } else {
      assert(0);
    }
  }

  void conducting_wall_J_lo(MfieldsState& mflds, int p, int d)
  {
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

  void conducting_wall_J_hi(MfieldsState& mflds, int p, int d)
  {
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

  void radiative_H_lo(MfieldsState& mflds, int p, int d)
  {
    set_lower_ghosts_to_nan(mflds, p, d, HX, false);

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
      if (radiation) {
        Real3 x3_s =
          (Real3(edge_idx) + Real3::unit(d1) * real_t(0.5)) * grid.domain.dx;
        Real3 x3_p =
          (Real3(edge_idx) + Real3::unit(d2) * real_t(0.5)) * grid.domain.dx;
        s = radiation->pulse_s_lower(grid.time(), d0, p, x3_s);
        p = radiation->pulse_p_lower(grid.time(), d0, p, x3_p);
      }

      F(H2, i3) =
        (4.f * s - 2.f * (F(E1, edge_idx) - background_e[d1]) -
         dtdx[d2] * (F(H0, edge_idx) - F(H0, edge_idx - Int3::unit(d2))) -
         (1.f - dtdx[d0]) * (F(H2, edge_idx) - background_h[d2]) +
         dt * F(J1, edge_idx)) /
          (1.f + dtdx[d0]) +
        background_h[d2];
      F(H1, i3) =
        (-4.f * p + 2.f * (F(E2, edge_idx) - background_e[d2]) -
         dtdx[d1] * (F(H0, edge_idx) - F(H0, edge_idx - Int3::unit(d1))) -
         (1.f - dtdx[d0]) * (F(H1, edge_idx) - background_h[d1]) +
         dt * F(J2, edge_idx)) /
          (1.f + dtdx[d0]) +
        background_h[d1];
    }
  }

  void radiative_H_hi(MfieldsState& mflds, int p, int d)
  {
    set_upper_ghosts_to_nan(mflds, p, d, HX, false);

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
      if (radiation) {
        Real3 x3_s =
          (Real3(edge_idx) + Real3::unit(d1) * real_t(0.5)) * grid.domain.dx;
        Real3 x3_p =
          (Real3(edge_idx) + Real3::unit(d2) * real_t(0.5)) * grid.domain.dx;
        s = radiation->pulse_s_upper(grid.time(), d0, p, x3_s);
        p = radiation->pulse_p_upper(grid.time(), d0, p, x3_p);
      }

      F(H2, i3) = (-4.f * s + 2.f * (F(E1, i3) - background_e[d1]) +
                   dtdx[d2] * (F(H0, i3) - F(H0, i3 - Int3::unit(d2))) -
                   (1.f - dtdx[d0]) * (F(H2, edge_idx) - background_h[d2]) -
                   dt * F(J1, i3)) /
                    (1.f + dtdx[d0]) +
                  background_h[d2];
      F(H1, i3) = (4.f * p - 2.f * (F(E2, i3) - background_e[d2]) +
                   dtdx[d1] * (F(H0, i3) - F(H0, i3 - Int3::unit(d1))) -
                   (1.f - dtdx[d0]) * (F(H1, edge_idx) - background_h[d1]) -
                   dt * F(J2, i3)) /
                    (1.f + dtdx[d0]) +
                  background_h[d1];
    }
  }

  Vec3<real_t> background_e = {0.0, 0.0, 0.0};
  Vec3<real_t> background_h = {0.0, 0.0, 0.0};

  RadiatingBoundary<real_t>* radiation;
};

// ======================================================================
// BndFieldsNone

// used by CUDA

template <class MFIELDS_STATE>
struct BndFieldsNone : BndFieldsBase
{
  using MfieldsState = MFIELDS_STATE;

  // clang-format off
  void fill_ghosts_E(MfieldsState& mflds) {};
  void fill_ghosts_H(MfieldsState& mflds) {};
  void add_ghosts_J(MfieldsState& mflds) {};
  // clang-format on
};
