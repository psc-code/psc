
#pragma once

#include <algorithm>

// ======================================================================
// SetupFields

template<typename MF>
struct SetupFields
{
  using Mfields = MF;
  using Fields = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;

  template<typename FUNC>
  static void set(Mfields& mf, FUNC func)
  {
    for (int p = 0; p < mf.n_patches(); ++p) {
      auto& patch = mf.grid().patches[p];
      Fields F(mf[p]);

      int n_ghosts = std::max({mf.ibn()[0], mf.ibn()[1], mf.ibn()[2]}); // FIXME, not pretty
      // FIXME, do we need the ghost points?
      mf.grid().Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
	  double x_nc = patch.x_nc(jx), y_nc = patch.y_nc(jy), z_nc = patch.z_nc(jz);
	  double x_cc = patch.x_cc(jx), y_cc = patch.y_cc(jy), z_cc = patch.z_cc(jz);
	  
	  double ncc[3] = { x_nc, y_cc, z_cc };
	  double cnc[3] = { x_cc, y_nc, z_cc };
	  double ccn[3] = { x_cc, y_cc, z_nc };
	  
	  double cnn[3] = { x_cc, y_nc, z_nc };
	  double ncn[3] = { x_nc, y_cc, z_nc };
	  double nnc[3] = { x_nc, y_nc, z_cc };
	  
	  F(HX, jx,jy,jz) += func(HX, ncc);
	  F(HY, jx,jy,jz) += func(HY, cnc);
	  F(HZ, jx,jy,jz) += func(HZ, ccn);
	  
	  F(EX, jx,jy,jz) += func(EX, cnn);
	  F(EY, jx,jy,jz) += func(EY, ncn);
	  F(EZ, jx,jy,jz) += func(EZ, nnc);
	  
	  F(JXI, jx,jy,jz) += func(JXI, cnn);
	  F(JYI, jx,jy,jz) += func(JYI, ncn);
	  F(JZI, jx,jy,jz) += func(JZI, nnc);
	});
    }
  }

  template<typename FUNC>
  static void set(MfieldsBase& mflds_base, FUNC func)
  {
    auto& mflds = mflds_base.get_as<Mfields>(0, 0);
    set(mflds, func);
    mflds_base.put_as(mflds, JXI, HX + 3);
  }
};
