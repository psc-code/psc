
#pragma once

#include "setup_fields.hxx"

// FIXME, duplicated only because vpic uses different comp numbers

#include <algorithm>

// ======================================================================
// SetupFields

using VpicConfig = VpicConfigPsc;

template <>
struct detail::SetupFields<VpicConfig::MfieldsState>
{
  using MfieldsState = VpicConfig::MfieldsState;

  template <typename FUNC>
  static void run(MfieldsState& mf, FUNC func)
  {
    for (int p = 0; p < mf.n_patches(); ++p) {
      auto& patch = mf.grid().patches[p];
      auto F = make_Fields3d<dim_xyz>(mf[p]);

      int n_ghosts =
        std::max({mf.ibn()[0], mf.ibn()[1], mf.ibn()[2]}); // FIXME, not pretty
      // FIXME, do we need the ghost points?
      mf.grid().Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
        double x_nc = patch.x_nc(jx), y_nc = patch.y_nc(jy),
               z_nc = patch.z_nc(jz);
        double x_cc = patch.x_cc(jx), y_cc = patch.y_cc(jy),
               z_cc = patch.z_cc(jz);

        double ncc[3] = {x_nc, y_cc, z_cc};
        double cnc[3] = {x_cc, y_nc, z_cc};
        double ccn[3] = {x_cc, y_cc, z_nc};

        double cnn[3] = {x_cc, y_nc, z_nc};
        double ncn[3] = {x_nc, y_cc, z_nc};
        double nnc[3] = {x_nc, y_nc, z_cc};

        F(MfieldsState::BX, jx, jy, jz) += func(HX, ncc);
        F(MfieldsState::BY, jx, jy, jz) += func(HY, cnc);
        F(MfieldsState::BZ, jx, jy, jz) += func(HZ, ccn);

        F(MfieldsState::EX, jx, jy, jz) += func(EX, cnn);
        F(MfieldsState::EY, jx, jy, jz) += func(EY, ncn);
        F(MfieldsState::EZ, jx, jy, jz) += func(EZ, nnc);
      });
    }
  }
};
