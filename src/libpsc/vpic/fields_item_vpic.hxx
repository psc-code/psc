
#pragma once

struct Item_vpic_fields
{
  using Mfields = MfieldsSingle;
  using Fields = Fields3d<typename Mfields::fields_t>;

  constexpr static const char* name = "fields_vpic";
  constexpr static int n_comps = 16;
  constexpr static fld_names_t fld_names()
  {
    return { "jx_ec", "jy_ec", "jz_ec",
	     "ex_ec", "ey_ec", "ez_ec",
	     "hx_fc", "hy_fc", "hz_fc",
	     "tcax_ec", "tcay_ec", "tcaz_ec",
	     "div_e_err_nc", "div_b_err_cc",
	     "rhob_nc", "rhof_nc", };
  }

  static void run(Mfields& mflds, Mfields& mres)
  {
    auto& grid = mflds.grid();
    
    for (int p = 0; p < mres.n_patches(); p++) {
      Fields F(mflds[p]), R(mres[p]);
      grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
	  for (int m = 0; m < 16; m++) {
	    R(m, ix,iy,iz) = F(m, ix,iy,iz);
	  }
	});
    }
  }
};

