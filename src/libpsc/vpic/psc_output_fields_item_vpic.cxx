
#include "psc_output_fields_item_private.h"

#include "vpic_config.h"
#include "fields.hxx"
#include "fields_item.hxx"
#include "psc_fields_as_single.h"
#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"

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
    for (int p = 0; p < mres.n_patches(); p++) {
      Fields F(mflds[p]), R(mres[p]);
      psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	for (int m = 0; m < 16; m++) {
	  R(m, ix,iy,iz) = F(m, ix,iy,iz);
	}
      } foreach_3d_end;
    }
  }
};

FieldsItemOps<FieldsItemFields<Item_vpic_fields>> psc_output_fields_item_vpic_fields_ops;

// ----------------------------------------------------------------------
// Moment_vpic_hydro

struct Moment_vpic_hydro : ItemMomentCRTP<Moment_vpic_hydro, MfieldsSingle>
{
  using Base = ItemMomentCRTP<Moment_vpic_hydro, MfieldsSingle>;
  using Mfields = MfieldsSingle;
  using Mparticles = MparticlesVpic;
  using Fields = Fields3d<typename Mfields::fields_t>;
  
  constexpr static char const* name = "hydro";
  constexpr static int n_comps = VPIC_HYDRO_N_COMP;
  constexpr static fld_names_t fld_names()
  {
    return { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
	     "px_nc", "py_nc", "pz_nc", "ke_nc",
	     "txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
	     "tzx_nc", "txy_nc", "_pad0", "_pad1", };
  }
  constexpr static int flags = POFI_BY_KIND;

  Moment_vpic_hydro(MPI_Comm comm)
    : Base(comm)
  {}
  
  void run(MparticlesVpic& mprts)
  {
    const auto& kinds = mprts.grid().kinds;
    auto& mres = *PscMfields<Mfields>{this->mres_}.sub();
    auto mf_hydro = MfieldsVpic{ppsc->grid(), 16, { 1, 1, 1 }};
    Simulation *sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    for (int kind = 0; kind < kinds.size(); kind++) {
      HydroArray *vmflds_hydro = mf_hydro.vmflds_hydro;
      Simulation_moments_run(sim, vmflds_hydro, mprts.vmprts, kind);
      
      auto& mf = mf_hydro.get_as<Mfields>(0, VPIC_HYDRO_N_COMP);
      for (int p = 0; p < mres.n_patches(); p++) {
	Fields F(mf[p]), R(mres[p]);
	psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	  for (int m = 0; m < VPIC_HYDRO_N_COMP; m++) {
	    R(m + kind * VPIC_HYDRO_N_COMP, ix,iy,iz) = F(m, ix,iy,iz);
	  }
	} foreach_3d_end;
      }
      mf_hydro.put_as(mf, 0, 0);
    }
  }
};

FieldsItemOps<FieldsItemMoment<Moment_vpic_hydro>> psc_output_fields_item_vpic_hydro_ops;

