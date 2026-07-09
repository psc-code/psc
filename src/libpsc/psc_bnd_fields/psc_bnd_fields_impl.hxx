#pragma once

#include "psc.h"
#include "kg/VecRange.hxx"
#include "fields.hxx"
#include "bnd_fields.hxx"
#include "radiating_bnd.hxx"
#include "field_bc_util.hxx"

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
              break;
            }
            case BND_FLD_OPEN: {
              psc::bnd::field::detail::set_lower_ghosts<dim_t>(
                mflds, p, d, EX, background_e_lo, false);
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
              break;
            }
            case BND_FLD_OPEN: {
              psc::bnd::field::detail::set_upper_ghosts<dim_t>(
                mflds, p, d, EX, background_e_hi, false);
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

  void conducting_wall_H_lo(MfieldsState& mflds, int p, int d)
  {
    psc::bnd::field::detail::set_lower_ghosts_to_nan<dim_t>(mflds, p, d, HX,
                                                            false);

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
    psc::bnd::field::detail::set_upper_ghosts_to_nan<dim_t>(mflds, p, d, HX,
                                                            false);

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
      if (radiation) {
        Real3 x3_s = (Real3(edge_idx) + Real3::unit(d1) * real_t(0.5)) *
                     Real3(grid.domain.dx);
        Real3 x3_p = (Real3(edge_idx) + Real3::unit(d2) * real_t(0.5)) *
                     Real3(grid.domain.dx);
        s = radiation->pulse_s_lower(grid.time(), d0, p, x3_s);
        p = radiation->pulse_p_lower(grid.time(), d0, p, x3_p);
      }

      F(H2, i3) =
        (4.f * s - 2.f * (F(E1, edge_idx) - background_e_lo[d1]) -
         dtdx[d2] * (F(H0, edge_idx) - F(H0, edge_idx - Int3::unit(d2))) -
         (1.f - dtdx[d0]) * (F(H2, edge_idx) - background_h_lo[d2]) +
         dt * F(J1, edge_idx)) /
          (1.f + dtdx[d0]) +
        background_h_lo[d2];
      F(H1, i3) =
        (-4.f * p + 2.f * (F(E2, edge_idx) - background_e_lo[d2]) -
         dtdx[d1] * (F(H0, edge_idx) - F(H0, edge_idx - Int3::unit(d1))) -
         (1.f - dtdx[d0]) * (F(H1, edge_idx) - background_h_lo[d1]) +
         dt * F(J2, edge_idx)) /
          (1.f + dtdx[d0]) +
        background_h_lo[d1];
    }
  }

  void radiative_H_hi(MfieldsState& mflds, int p, int d)
  {
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
      if (radiation) {
        Real3 x3_s = (Real3(edge_idx) + Real3::unit(d1) * real_t(0.5)) *
                     Real3(grid.domain.dx);
        Real3 x3_p = (Real3(edge_idx) + Real3::unit(d2) * real_t(0.5)) *
                     Real3(grid.domain.dx);
        s = radiation->pulse_s_upper(grid.time(), d0, p, x3_s);
        p = radiation->pulse_p_upper(grid.time(), d0, p, x3_p);
      }

      F(H2, i3) = (-4.f * s + 2.f * (F(E1, i3) - background_e_hi[d1]) +
                   dtdx[d2] * (F(H0, i3) - F(H0, i3 - Int3::unit(d2))) -
                   (1.f - dtdx[d0]) * (F(H2, edge_idx) - background_h_hi[d2]) -
                   dt * F(J1, i3)) /
                    (1.f + dtdx[d0]) +
                  background_h_hi[d2];
      F(H1, i3) = (4.f * p - 2.f * (F(E2, i3) - background_e_hi[d2]) +
                   dtdx[d1] * (F(H0, i3) - F(H0, i3 - Int3::unit(d1))) -
                   (1.f - dtdx[d0]) * (F(H1, edge_idx) - background_h_hi[d1]) -
                   dt * F(J2, i3)) /
                    (1.f + dtdx[d0]) +
                  background_h_hi[d1];
    }
  }

  Vec3<real_t> background_e_lo = {0.0, 0.0, 0.0};
  Vec3<real_t> background_e_hi = {0.0, 0.0, 0.0};

  Vec3<real_t> background_h_lo = {0.0, 0.0, 0.0};
  Vec3<real_t> background_h_hi = {0.0, 0.0, 0.0};

  RadiatingBoundary<real_t>* radiation = nullptr;
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
