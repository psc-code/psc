
#include "psc_output_fields_item_private.h"

#include "vpic_config.h"
#include "fields.hxx"
#include "fields_item.hxx"
#include "psc_fields_as_c.h"
#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"

using fields_t = mfields_t::fields_t;
using Fields = Fields3d<fields_t>;

struct Item_vpic_fields
{
  using mfields_t = mfields_t;
  constexpr static const char* name = "vpic_fields";
  constexpr static int n_comps = 16;
  constexpr static fld_names_t fld_names()
  {
#if 0
    return { "ex_ec", "ey_ec", "ez_ec", "dive_nc",
	     "hx_fc", "hy_fc", "hz_fc", "divb_cc",
	     "tcax", "tcay", "tcaz", "rhob_nc",
	     "jx_ec", "jy_ec", "jz_ec", "rho_nc", };
#else
    return { "jx_ec", "jy_ec", "jz_ec",
	     "ex_ec", "ey_ec", "ez_ec",
	     "hx_fc", "hy_fc", "hz_fc",
	     "tcax_ec", "tcay_ec", "tcaz_ec",
	     "div_e_err_nc", "div_b_err_cc",
	     "rhob_nc", "rhof_nc", };
#endif
  }

  static void run(mfields_t mflds, mfields_t mres)
  {
    for (int p = 0; p < mres->n_patches(); p++) {
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

struct Moment_vpic_hydro
{
  using mfields_t = mfields_t;
  using mparticles_t = PscMparticlesVpic;
  
  constexpr static char const* name = "vpic_hydro";
  constexpr static int n_comps = VPIC_HYDRO_N_COMP;
  constexpr static fld_names_t fld_names()
  {
    return { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
	     "px_nc", "py_nc", "pz_nc", "ke_nc",
	     "txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
	     "tzx_nc", "txy_nc", "_pad0", "_pad1", };
  }
  constexpr static int flags = POFI_BY_KIND;
  
  static void run(mfields_t mres, mparticles_t mprts)
  {
    struct psc_mfields *mflds_hydro = psc_mfields_create(psc_mfields_comm(mres.mflds()));
    psc_mfields_set_type(mflds_hydro, "vpic");
    psc_mfields_set_param_int(mflds_hydro, "nr_fields", 16);
    int ibn[3] = { 1, 1, 1 };
    psc_mfields_set_param_int3(mflds_hydro, "ibn", ibn);
    mflds_hydro->grid = &ppsc->grid();
    psc_mfields_setup(mflds_hydro);
    PscMfieldsVpic mf_hydro(mflds_hydro);
    
    Simulation *sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    for (int kind = 0; kind < ppsc->nr_kinds; kind++) {
      HydroArray *vmflds_hydro = mf_hydro->vmflds_hydro;
      Simulation_moments_run(sim, vmflds_hydro, mprts->vmprts, kind);
      
      mfields_t mf = mflds_hydro->get_as<mfields_t>(0, VPIC_HYDRO_N_COMP);
      for (int p = 0; p < mres->n_patches(); p++) {
	Fields F(mf[p]), R(mres[p]);
	psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	  for (int m = 0; m < VPIC_HYDRO_N_COMP; m++) {
	    R(m + kind * VPIC_HYDRO_N_COMP, ix,iy,iz) = F(m, ix,iy,iz);
	  }
	} foreach_3d_end;
      }
      mf.put_as(mflds_hydro, 0, 0);
    }
    
    psc_mfields_destroy(mflds_hydro);
  }
};

FieldsItemOps<ItemMoment<Moment_vpic_hydro>> psc_output_fields_item_vpic_hydro_ops;

