
#pragma once

// ----------------------------------------------------------------------
// OutputFieldsVpic

struct OutputFieldsVpic
{
  struct Result
  {
    MfieldsStateVpic& mflds;
    const char* name;
    std::vector<std::string> comp_names;
  };
  
  Result operator()(MfieldsStateVpic& mflds)
  {
    std::vector<std::string> comp_names = { "ex_ec", "ey_ec", "ez_ec", "div_e_err_nc",
					    "hx_fc", "hy_fc", "hz_fc", "div_b_err_cc"
					    "tcax_ec", "tcay_ec", "tcaz_ec", "rhob_nc",
					    "jx_ec", "jy_ec", "jz_c", "rhof_nc" };
    return {mflds, "fields_vpic", comp_names};
  }
};

// ----------------------------------------------------------------------
// OutputHydroVpic

struct OutputHydroVpic
{
  struct Result
  {
    MfieldsSingle& mflds;
    const char* name;
    std::vector<std::string> comp_names;
  };
  
  constexpr static fld_names_t fld_names() {
    return { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
	"px_nc", "py_nc", "pz_nc", "ke_nc",
	"txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
	"tzx_nc", "txy_nc", "_pad0", "_pad1" };
  }

  OutputHydroVpic(const Grid_t& grid)
    : mflds_res_{grid, MfieldsHydroVpic::N_COMP * int(grid.kinds.size()), {1,1,1}}
  {}

  Result operator()(MparticlesVpic& mprts, MfieldsHydroVpic& mflds_hydro, Interpolator& interpolator)
  {
    // This relies on load_interpolator_array() having been called earlier

    auto& grid = mflds_res_.grid();
    Simulation *sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    Particles& vmprts = mprts.vmprts_;
    HydroArray& vmflds = *mflds_hydro.vmflds_hydro;

    std::vector<std::string> comp_names;
    comp_names.reserve(mflds_res_.n_comps());

    
    for (int kind = 0; kind < grid.kinds.size(); ++kind) {
      vmflds.clear();
      
      // FIXME, just iterate over species instead?
      typename Particles::const_iterator sp = vmprts.find(kind);
      vmprts.accumulate_hydro_p(vmflds, sp, interpolator);
      
      vmflds.synchronize();
      
      for (int p = 0; p < mflds_res_.n_patches(); p++) {
	auto res = mflds_res_[p];
	auto H = mflds_hydro[p];
	grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	    for (int m = 0; m < MfieldsHydroVpic::N_COMP; m++) {
	      res(m + kind * MfieldsHydroVpic::N_COMP, i,j,k) = H(m, i,j,k);
	    }
	  });
      }

      for (int m = 0; m < MfieldsHydroVpic::N_COMP; m++) {
	comp_names.emplace_back(std::string(fld_names()[m]) + "_" + grid.kinds[kind].name);
      }
    }

    return {mflds_res_, "hydro_vpic", comp_names};
  }

private:
  MfieldsSingle mflds_res_;
};

