
#include "psc.h"
#include "kg/VecRange.hxx"
#include "fields.hxx"
#include "bnd_fields.hxx"

#include <mrc_bits.h>

#include <limits>

// #define DEBUG

template <typename MFIELDS_STATE, typename Dim>
struct BndFields_ : BndFieldsBase
{
  using Self = BndFields_<MFIELDS_STATE, Dim>;
  using MfieldsState = MFIELDS_STATE;
  using real_t = typename MfieldsState::real_t;
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
              set_background_E_lo(mflds, p, d);
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
              set_background_E_hi(mflds, p, d);
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
              add_background_H_lo(mflds, p, d);
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
              conducting_wall_H_hi(mflds, p, d);
              break;
            }
            case BND_FLD_OPEN: {
              radiative_H_hi(mflds, p, d);
              add_background_H_hi(mflds, p, d);
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

  static void fields_t_set_nan(real_t* f)
  {
    *f = std::numeric_limits<real_t>::quiet_NaN();
  }

  void conducting_wall_E_lo(MfieldsState& mflds, int p, int d)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
#ifdef DEBUG
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, F.ib_[0]);
             ix < std::min(ldims[0] + 2, F.ib_[0] + F.im_[0]); ix++) {
          fields_t_set_nan(&F(EX, ix, -1, iz));
          fields_t_set_nan(&F(EX, ix, -2, iz));
          fields_t_set_nan(&F(EY, ix, -1, iz));
          fields_t_set_nan(&F(EY, ix, -2, iz));
          fields_t_set_nan(&F(EZ, ix, -1, iz));
          fields_t_set_nan(&F(EZ, ix, -2, iz));
        }
      }
#endif
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
#ifdef DEBUG
      for (int iy = -2; iy < ldims[1] + 2; iy++) {
        for (int ix = -2; ix < ldims[0] + 2; ix++) {
          fields_t_set_nan(&F(EX, ix, iy, -1));
          fields_t_set_nan(&F(EX, ix, iy, -2));
          fields_t_set_nan(&F(EY, ix, iy, -1));
          fields_t_set_nan(&F(EY, ix, iy, -2));
          fields_t_set_nan(&F(EZ, ix, iy, -1));
          fields_t_set_nan(&F(EZ, ix, iy, -2));
        }
      }
#endif
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
    auto F = make_Fields3d<dim_t>(mflds[p]);
    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
      int my _mrc_unused = ldims[1];
#ifdef DEBUG
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, F.ib_[0]);
             ix < std::min(ldims[0] + 2, F.ib_[0] + F.im_[0]); ix++) {
          fields_t_set_nan(&F(EX, ix, my, iz));
          fields_t_set_nan(&F(EX, ix, my + 1, iz));
          fields_t_set_nan(&F(EY, ix, my, iz));
          fields_t_set_nan(&F(EY, ix, my + 1, iz));
          fields_t_set_nan(&F(EZ, ix, my, iz));
          fields_t_set_nan(&F(EZ, ix, my + 1, iz));
        }
      }
#endif
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
#ifdef DEBUG
      for (int iy = -2; iy < ldims[1] + 2; iy++) {
        for (int ix = -2; ix < ldims[0] + 2; ix++) {
          fields_t_set_nan(&F(EX, ix, iy, mz));
          fields_t_set_nan(&F(EX, ix, iy, mz + 1));
          fields_t_set_nan(&F(EY, ix, iy, mz));
          fields_t_set_nan(&F(EY, ix, iy, mz + 1));
          fields_t_set_nan(&F(EZ, ix, iy, mz));
          fields_t_set_nan(&F(EZ, ix, iy, mz + 1));
        }
      }
#endif
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
    auto F = make_Fields3d<dim_t>(mflds[p]);
    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
#ifdef DEBUG
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, F.ib_[0]);
             ix < std::min(ldims[0] + 2, F.ib_[0] + F.im_[0]); ix++) {
          fields_t_set_nan(&F(HX, ix, -1, iz));
          fields_t_set_nan(&F(HX, ix, -2, iz));
          fields_t_set_nan(&F(HY, ix, -1, iz));
          fields_t_set_nan(&F(HY, ix, -2, iz));
          fields_t_set_nan(&F(HZ, ix, -1, iz));
          fields_t_set_nan(&F(HZ, ix, -2, iz));
        }
      }
#endif
      for (int iz = -1; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, -1, iz) = -F(HX, ix, 0, iz);

          F(HY, ix, -1, iz) = F(HY, ix, 1, iz);

          F(HZ, ix, -1, iz) = -F(HZ, ix, 0, iz);
        }
      }
    } else if (d == 2) {
#ifdef DEBUG
      for (int iy = -2; iy < ldims[1] + 2; iy++) {
        for (int ix = -2; ix < ldims[0] + 2; ix++) {
          fields_t_set_nan(&F(HX, ix, iy, -1));
          fields_t_set_nan(&F(HX, ix, iy, -2));
          fields_t_set_nan(&F(HY, ix, iy, -1));
          fields_t_set_nan(&F(HY, ix, iy, -2));
          fields_t_set_nan(&F(HZ, ix, iy, -1));
          fields_t_set_nan(&F(HZ, ix, iy, -2));
        }
      }
#endif
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
    auto F = make_Fields3d<dim_t>(mflds[p]);

    const int* ldims = mflds.grid().ldims;
    Int3 ib = mflds.ib(), im = mflds.im();

    if (d == 1) {
      int my _mrc_unused = ldims[1];
#ifdef DEBUG
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, F.ib_[0]);
             ix < std::min(ldims[0] + 2, F.ib_[0] + F.im_[0]); ix++) {
          fields_t_set_nan(&F(HX, ix, my, iz));
          fields_t_set_nan(&F(HX, ix, my + 1, iz));
          fields_t_set_nan(&F(HY, ix, my + 1, iz));
          fields_t_set_nan(&F(HZ, ix, my, iz));
          fields_t_set_nan(&F(HZ, ix, my + 1, iz));
        }
      }
#endif
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
#ifdef DEBUG
      for (int iy = -2; iy < ldims[1] + 2; iy++) {
        for (int ix = -2; ix < ldims[0] + 2; ix++) {
          fields_t_set_nan(&F(HX, ix, iy, mz));
          fields_t_set_nan(&F(HX, ix, iy, mz + 1));
          fields_t_set_nan(&F(HY, ix, iy, mz));
          fields_t_set_nan(&F(HY, ix, iy, mz + 1));
          fields_t_set_nan(&F(HZ, ix, iy, mz + 1));
        }
      }
#endif
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
    auto F = make_Fields3d<dim_t>(mflds[p]);
    const Grid_t& grid = mflds.grid();
    Int3 ldims = grid.ldims;
    Int3 ib = mflds.ib();
    Int3 im = mflds.im();
    real_t dt = grid.dt;
    real_t dx = grid.domain.dx[0];
    real_t dy = grid.domain.dx[1];
    real_t dz = grid.domain.dx[2];

    if (d == 1) {
#ifdef DEBUG
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          fields_t_set_nan(&F(HX, ix, -1, iz));
          fields_t_set_nan(&F(HX, ix, -2, iz));
          fields_t_set_nan(&F(HY, ix, -1, iz));
          fields_t_set_nan(&F(HY, ix, -2, iz));
          fields_t_set_nan(&F(HZ, ix, -1, iz));
          fields_t_set_nan(&F(HZ, ix, -2, iz));
        }
      }
#endif
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, -1, iz) =
            (/* + 4.f * C_s_pulse_y1(x,y,z+0.5*dz,t) */
             -2.f * F(EZ, ix, 0, iz)
             /*- dt/dx * (F(HY, ix,0,iz) - F(HY, ix-1,0,iz)) */
             - (1.f - dt / dy) * F(HX, ix, 0, iz) + dt * F(JZI, ix, 0, iz)) /
            (1.f + dt / dy);
          F(HZ, ix, -1, iz) =
            (/* + 4.f * C_p_pulse_y1(x+.5*dx,y,z,t) */
             +2.f * F(EX, ix, 0, iz) -
             dt / dz * (F(HY, ix, 0, iz) - F(HY, ix, 0, iz - 1)) -
             (1.f - dt / dy) * F(HZ, ix, 0, iz) + dt * F(JXI, ix, 0, iz)) /
            (1.f + dt / dy);
        }
      }
    } else {
      assert(0);
    }
  }

  void radiative_H_hi(MfieldsState& mflds, int p, int d)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    const Grid_t& grid = mflds.grid();
    Int3 ldims = grid.ldims;
    Int3 ib = mflds.ib();
    Int3 im = mflds.im();
    real_t dt = grid.dt;
    real_t dx = grid.domain.dx[0];
    real_t dy = grid.domain.dx[1];
    real_t dz = grid.domain.dx[2];

    if (d == 1) {
      int my _mrc_unused = ldims[1];
#ifdef DEBUG
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          fields_t_set_nan(&F(HX, ix, my, iz));
          fields_t_set_nan(&F(HX, ix, my + 1, iz));
          fields_t_set_nan(&F(HY, ix, my, iz));
          fields_t_set_nan(&F(HY, ix, my + 1, iz));
          fields_t_set_nan(&F(HZ, ix, my + 1, iz));
        }
      }
#endif
      for (int iz = -2; iz < ldims[2] + 2; iz++) {
        for (int ix = std::max(-2, ib[0]);
             ix < std::min(ldims[0] + 2, ib[0] + im[0]); ix++) {
          F(HX, ix, my, iz) =
            (/* + 4.f * C_s_pulse_y2(x,y,z+0.5*dz,t) */
             +2.f * F(EZ, ix, my, iz)
             /*+ dt/dx * (F(HY, ix,my,iz) - F(HY, ix-1,my,iz)) */
             - (1.f - dt / dy) * F(HX, ix, my - 1, iz) -
             dt * F(JZI, ix, my, iz)) /
            (1.f + dt / dy);
          F(HZ, ix, my, iz) =
            (/* + 4.f * C_p_pulse_y2(x+.5*dx,y,z,t) */
             -2.f * F(EX, ix, my, iz) +
             dt / dz * (F(HY, ix, my, iz) - F(HY, ix, my, iz - 1)) -
             (1.f - dt / dy) * F(HZ, ix, my - 1, iz) -
             dt * F(JXI, ix, my, iz)) /
            (1.f + dt / dy);
        }
      }
    } else {
      assert(0);
    }
  }

  void set_background_E_lo(MfieldsState& mflds, int p, int d)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    Int3 start = mflds.ib();
    Int3 stop = mflds.im();
    stop[d] = 0;
    for (Int3 i3 : VecRange(start, stop)) {
      F(EX, i3) = background_e[0];
      F(EY, i3) = background_e[1];
      F(EZ, i3) = background_e[2];
    }
  }

  void set_background_E_hi(MfieldsState& mflds, int p, int d)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    Int3 start = mflds.ib();
    Int3 stop = mflds.im();
    start[d] = mflds.grid().ldims[d] + 1;

    Int3 neg1 = {0, 0, 0};
    neg1[d] = -1;

    for (Int3 i3 : VecRange(start, stop)) {
      F(EX, i3) = background_e[0];
      F(EY, i3) = background_e[1];
      F(EZ, i3) = background_e[2];

      // the other two components of e are in the domain at this index
      F(EX + d, i3 + neg1) = background_e[d];
    }
  }

  void add_background_H_lo(MfieldsState& mflds, int p, int d)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    Int3 start = mflds.ib();
    Int3 stop = mflds.im();
    stop[d] = 0;
    for (Int3 i3 : VecRange(start, stop)) {
      F(HX, i3) += background_h[0];
      F(HY, i3) += background_h[1];
      F(HZ, i3) += background_h[2];
    }
  }

  void add_background_H_hi(MfieldsState& mflds, int p, int d)
  {
    auto F = make_Fields3d<dim_t>(mflds[p]);
    Int3 start = mflds.ib();
    Int3 stop = mflds.im();
    start[d] = mflds.grid().ldims[d] + 1;

    Int3 neg1 = {0, 0, 0};
    neg1[d] = -1;

    for (Int3 i3 : VecRange(start, stop)) {
      F(HX, i3) += background_h[0];
      F(HY, i3) += background_h[1];
      F(HZ, i3) += background_h[2];

      // the third component of h is in the domain at this index
      int d1 = (d + 1) % 3;
      int d2 = (d + 2) % 3;
      F(HX + d1, i3 + neg1) += background_h[d1];
      F(HX + d2, i3 + neg1) += background_h[d2];
    }
  }

  Vec3<real_t> background_e = {0.0, 0.0, 0.0};
  Vec3<real_t> background_h = {0.0, 0.0, 0.0};
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
