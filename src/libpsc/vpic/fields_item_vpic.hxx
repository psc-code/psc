
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

struct Moment_vpic_hydro : ItemMomentCRTP<Moment_vpic_hydro, MfieldsSingle>
{
  using Base = ItemMomentCRTP<Moment_vpic_hydro, MfieldsSingle>;
  using Mfields = MfieldsSingle;
  using MfieldsState = MfieldsStateVpic;
  using Mparticles = MparticlesVpic;
  using Fields = Fields3d<typename Mfields::fields_t>;
  
  constexpr static char const* name = "hydro";
  constexpr static int n_comps = MfieldsHydroVpic::N_COMP;
  constexpr static fld_names_t fld_names()
  {
    return { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
	     "px_nc", "py_nc", "pz_nc", "ke_nc",
	     "txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
	     "tzx_nc", "txy_nc", "_pad0", "_pad1", };
  }
  constexpr static int flags = POFI_BY_KIND;

  Moment_vpic_hydro(const Grid_t& grid, MPI_Comm comm)
    : Base(grid, comm)
  {}
  
  void run(MparticlesVpic& mprts)
  {
    const auto& kinds = mprts.grid().kinds;
    auto& mres = this->mres_;
    auto& grid = mprts.grid();
    auto mf_hydro = MfieldsHydroVpic{ppsc->grid(), 16, { 1, 1, 1 }};
    Simulation *sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
      // This relies on load_interpolator_array() having been called earlier
    for (int kind = 0; kind < kinds.size(); kind++) {
      HydroArray& hydro = *mf_hydro.vmflds_hydro;
      Particles& vmprts = mprts.vmprts_;

      hydro.clear();

      // FIXME, just iterate over species instead?
      typename Particles::const_iterator sp = vmprts.find(kind);
      vmprts.accumulate_hydro_p(hydro, sp, *sim->interpolator_);
      
      hydro.synchronize();
      
      for (int p = 0; p < mres.n_patches(); p++) {
	Fields R(mres[p]);
	auto F = mf_hydro[p];
	grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
	    for (int m = 0; m < MfieldsHydroVpic::N_COMP; m++) {
	      R(m + kind * MfieldsHydroVpic::N_COMP, ix,iy,iz) = F(m, ix,iy,iz);
	    }
	  });
      }
    }
  }
};

