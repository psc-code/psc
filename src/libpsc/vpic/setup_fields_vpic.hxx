
#pragma once

// FIXME, duplicated only because vpic uses different comp numbers

#include <algorithm>

// ======================================================================
// SetupFields

template<>
struct SetupFields<MfieldsStateVpic>
{
  using Mfields = MfieldsStateVpic;
  using Fields = typename Mfields::fields_t;

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
	  
	  F(VPIC_MFIELDS_BX, jx,jy,jz) += func(HX, ncc);
	  F(VPIC_MFIELDS_BY, jx,jy,jz) += func(HY, cnc);
	  F(VPIC_MFIELDS_BZ, jx,jy,jz) += func(HZ, ccn);
	  
	  F(VPIC_MFIELDS_EX, jx,jy,jz) += func(EX, cnn);
	  F(VPIC_MFIELDS_EY, jx,jy,jz) += func(EY, ncn);
	  F(VPIC_MFIELDS_EZ, jx,jy,jz) += func(EZ, nnc);
	});
    }
  }
};
