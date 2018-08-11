
#pragma once

// ----------------------------------------------------------------------
// Item_vpic_fields

struct Item_vpic_fields
{
  using Mfields = MfieldsSingle;
  using MfieldsState = MfieldsStateVpic;
  using Fields = Fields3d<typename Mfields::fields_t>;
  using FieldsState = Fields3d<typename MfieldsState::fields_t>;

  constexpr static const char* name = "fields_vpic";
  constexpr static int n_comps = 16;
  constexpr static fld_names_t fld_names()
  {
    return { "ex_ec", "ey_ec", "ez_ec", "div_e_err_nc",
    	     "hx_fc", "hy_fc", "hz_fc", "div_b_err_cc",
	     "tcax_ec", "tcay_ec", "tcaz_ec", "rhob_nc",
	     "jx_ec", "jy_ec", "jz_c", "rhof_nc" };
  }

  static void run(MfieldsState& mflds, Mfields& mres)
  {
    auto& grid = mflds.grid();
    
    for (int p = 0; p < mres.n_patches(); p++) {
      FieldsState F(mflds[p]);
      Fields R(mres[p]);
      grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
	  for (int m = 0; m < 16; m++) {
	    R(m, ix,iy,iz) = F(m, ix,iy,iz);
	  }
	});
    }
  }
};

// ----------------------------------------------------------------------
// Moment_vpic_hydro

struct Moment_vpic_hydro
{
  using Mfields = MfieldsSingle;
  using Mparticles = MparticlesVpic;
  using Fields = Fields3d<typename Mfields::fields_t>;
  
  constexpr static fld_names_t fld_names() {
    return { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
	"px_nc", "py_nc", "pz_nc", "ke_nc",
	"txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
	"tzx_nc", "txy_nc", "_pad0", "_pad1" };
  }

  Moment_vpic_hydro(const Grid_t& grid)
    : mflds_res_{grid, MfieldsHydroVpic::N_COMP * int(grid.kinds.size()), {1,1,1}}
  {}

  std::vector<std::string> comp_names()
  {
    auto& grid = mflds_res_.grid();
    
    std::vector<std::string> comp_names;
    comp_names.reserve(mflds_res_.n_comps());

    for (int kind = 0; kind < grid.kinds.size(); kind++) {
      for (int m = 0; m < MfieldsHydroVpic::N_COMP; m++) {
	comp_names.emplace_back(std::string(fld_names()[m]) + "_" + grid.kinds[kind].name);
      }
    }
    return comp_names;
  }
  
  MfieldsSingle& run(MparticlesVpic& mprts, MfieldsHydroVpic& mflds_hydro)
  {
    // This relies on load_interpolator_array() having been called earlier

    auto& grid = mflds_res_.grid();
    Simulation *sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    Particles& vmprts = mprts.vmprts_;
    
    for (int kind = 0; kind < grid.kinds.size(); ++kind) {
      mflds_hydro.vmflds_hydro->clear();
      
      // FIXME, just iterate over species instead?
      typename Particles::const_iterator sp = vmprts.find(kind);
      vmprts.accumulate_hydro_p(*mflds_hydro.vmflds_hydro, sp, *sim->interpolator_);
      
      mflds_hydro.vmflds_hydro->synchronize();
      
      for (int p = 0; p < mflds_res_.n_patches(); p++) {
	auto res = mflds_res_[p];
	auto H = mflds_hydro[p];
	grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	    for (int m = 0; m < MfieldsHydroVpic::N_COMP; m++) {
	      res(m + kind * MfieldsHydroVpic::N_COMP, i,j,k) = H(m, i,j,k);
	    }
	  });
      }
    }

    return mflds_res_;
  }

private:
  MfieldsSingle mflds_res_;
};

