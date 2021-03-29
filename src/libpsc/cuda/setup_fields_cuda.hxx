
#pragma once

#include "setup_fields.hxx"
#include "psc_fields_cuda.h"

template <>
template <typename FUNC>
void detail::SetupFields<MfieldsStateCuda>::run(Mfields& mf, FUNC&& func)
{
  auto h_mf = hostMirror(mf);
  copy(mf, h_mf);
  for (int p = 0; p < mf.n_patches(); ++p) {
    auto& patch = mf.grid().patches[p];
    auto flds = make_Fields3d<dim_xyz>(h_mf[p]);

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

      flds(HX, jx, jy, jz) += func(HX, ncc);
      flds(HY, jx, jy, jz) += func(HY, cnc);
      flds(HZ, jx, jy, jz) += func(HZ, ccn);

      flds(EX, jx, jy, jz) += func(EX, cnn);
      flds(EY, jx, jy, jz) += func(EY, ncn);
      flds(EZ, jx, jy, jz) += func(EZ, nnc);

      flds(JXI, jx, jy, jz) += func(JXI, cnn);
      flds(JYI, jx, jy, jz) += func(JYI, ncn);
      flds(JZI, jx, jy, jz) += func(JZI, nnc);
    });
  }
  copy(h_mf, mf);
}
