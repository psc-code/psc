
#pragma once

#include "fields_item.hxx"
#include <psc_fields_single.h>
#ifdef USE_VPIC
#include "vpic_hydro_ops.hxx"
#endif
#include "psc_hydro_ops.hxx"

// ----------------------------------------------------------------------
// OutputFieldsVpic

template<typename MfieldsState>
struct OutputFieldsVpic
{
  struct Result
  {
    MfieldsState& mflds;
    const char* name;
    std::vector<std::string> comp_names;
  };
  
  Result operator()(MfieldsState& mflds)
  {
    std::vector<std::string> comp_names = { "ex_ec", "ey_ec", "ez_ec", "div_e_err_nc",
					    "hx_fc", "hy_fc", "hz_fc", "div_b_err_cc",
					    "tcax_ec", "tcay_ec", "tcaz_ec", "rhob_nc",
					    "jx_ec", "jy_ec", "jz_ec", "rhof_nc" };
    return {mflds, "fields_vpic", comp_names};
  }
};

// ----------------------------------------------------------------------
// OutputHydroVpic

template<typename HydroOps>
struct OutputHydroVpic_
{
  using Mparticles = typename HydroOps::Mparticles;
  using MfieldsHydro = typename HydroOps::MfieldsHydro;
  using MfieldsInterpolator = typename HydroOps::MfieldsInterpolator;
  
  struct Result
  {
    MfieldsSingle& mflds;
    const char* name;
    std::vector<std::string> comp_names;
  };
  
  static std::vector<std::string> fld_names() {
    return { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
	"px_nc", "py_nc", "pz_nc", "ke_nc",
	"txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
	"tzx_nc", "txy_nc", "_pad0", "_pad1" };
  }

  OutputHydroVpic_(const Grid_t& grid)
    : mflds_res_{grid, MfieldsHydro::N_COMP * int(grid.kinds.size()), {1,1,1}},
      kinds_{grid.kinds}
  {}

  Result operator()(Mparticles& mprts, MfieldsHydro& mflds_hydro, MfieldsInterpolator& interpolator)
  {
    // This relies on load_interpolator_array() having been called earlier

    std::vector<std::string> comp_names;
    comp_names.reserve(mflds_res_.n_comps());
    
    for (int kind = 0; kind < kinds_.size(); ++kind) {
      HydroOps::clear(mflds_hydro);
      
      auto prts = mprts[0];
      // FIXME, just iterate over species instead?
      auto sp = std::find_if(prts.cbegin(), prts.cend(),
			     [&](const typename Mparticles::Species& sp) { return sp.id == kind; });
      HydroOps::accumulate_hydro_p(mflds_hydro, sp, interpolator);
      
      HydroOps::synchronize(mflds_hydro);
      
      for (int p = 0; p < mflds_res_.n_patches(); p++) {
	auto res = mflds_res_[p];
	auto H = mflds_hydro[p];
	mflds_res_.Foreach_3d(0, 0, [&](int i, int j, int k) {
	    for (int m = 0; m < MfieldsHydro::N_COMP; m++) {
	      res(m + kind * MfieldsHydro::N_COMP, i,j,k) = H(m, i,j,k);
	    }
	  });
      }

      for (int m = 0; m < MfieldsHydro::N_COMP; m++) {
	comp_names.emplace_back(fld_names()[m] + "_" + kinds_[kind].name);
      }
    }

    return {mflds_res_, "hydro_vpic", comp_names};
  }

private:
  MfieldsSingle mflds_res_;
  Grid_t::Kinds kinds_;
};

#ifdef USE_VPIC
template<typename Mparticles, typename MfieldsHydro, typename MfieldsInterpolator>
using OutputHydroVpicWrap = OutputHydroVpic_<VpicHydroOps<Mparticles, MfieldsHydro, MfieldsInterpolator>>;
#endif

template<typename Mparticles, typename MfieldsHydro, typename MfieldsInterpolator>
using OutputHydroVpic = OutputHydroVpic_<PscHydroOps<Mparticles, MfieldsHydro, MfieldsInterpolator>>;
