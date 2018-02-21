
#include "psc_output_fields_item_private.h"

#include "vpic_config.h"
#include "fields.hxx"
#include "psc_fields_as_c.h"
#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"

using Fields = Fields3d<fields_t>;

// ----------------------------------------------------------------------
// run_all_vpic_fields

static void
run_all_vpic_fields(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		   struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  mfields_t mf = mflds_base->get_as<mfields_t>(0, 16);
  mfields_t mf_res(mres);
  
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      for (int m = 0; m < 16; m++) {
	R(m, ix,iy,iz) = F(m, ix,iy,iz);
      }
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_output_fields_item "vpic_fields"

struct psc_output_fields_item_ops_vpic : psc_output_fields_item_ops {
  psc_output_fields_item_ops_vpic() {
    name      = "vpic_fields";
    nr_comp   = 16;
#if 0
#if 0
    fld_names = { "ex_ec", "ey_ec", "ez_ec", "dive_nc",
		   "hx_fc", "hy_fc", "hz_fc", "divb_cc",
		   "tcax", "tcay", "tcaz", "rhob_nc",
		   "jx_ec", "jy_ec", "jz_ec", "rho_nc", };
#else
    fld_names = { "jx_ec", "jy_ec", "jz_ec",
		   "ex_ec", "ey_ec", "ez_ec",
		   "hx_fc", "hy_fc", "hz_fc",
		   "tcax_ec", "tcay_ec", "tcaz_ec",
		   "div_e_err_nc", "div_b_err_cc",
		   "rhob_nc", "rhof_nc", };
#endif
#endif
    run_all   = run_all_vpic_fields;
  }
} psc_output_fields_item_vpic_fields_ops;

// ----------------------------------------------------------------------
// run_all_vpic_hydro

static void
run_all_vpic_hydro(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		   struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  struct psc_mfields *mflds_hydro = psc_mfields_create(psc_mfields_comm(mres));
  psc_mfields_set_type(mflds_hydro, "vpic");
  psc_mfields_set_param_obj(mflds_hydro, "domain", mres->domain);
  psc_mfields_set_param_int(mflds_hydro, "nr_fields", 16);
  int ibn[3] = { 1, 1, 1 };
  psc_mfields_set_param_int3(mflds_hydro, "ibn", ibn);
  psc_mfields_setup(mflds_hydro);

  mparticles_vpic_t mprts = mprts_base->get_as<mparticles_vpic_t>();

  Simulation *sim;
  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);

  for (int kind = 0; kind < ppsc->nr_kinds; kind++) {
    HydroArray *vmflds_hydro = psc_mfields_vpic(mflds_hydro)->vmflds_hydro;
    Simulation_moments_run(sim, vmflds_hydro, mprts->vmprts, kind);
    
    mfields_t mf = mflds_hydro->get_as<mfields_t>(0, VPIC_HYDRO_N_COMP);
    mfields_t mf_res(mres);
    for (int p = 0; p < mf_res->n_patches(); p++) {
      Fields F(mf[p]), R(mf_res[p]);
      psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	for (int m = 0; m < VPIC_HYDRO_N_COMP; m++) {
	  R(m + kind * VPIC_HYDRO_N_COMP, ix,iy,iz) = F(m, ix,iy,iz);
	}
      } foreach_3d_end;
    }
    mf.put_as(mflds_hydro, 0, 0);
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);
  
  psc_mfields_destroy(mflds_hydro);
}

// ----------------------------------------------------------------------
// psc_output_fields_item "vpic_hydro"

struct psc_output_fields_item_ops_vpic_hydro : psc_output_fields_item_ops {
  psc_output_fields_item_ops_vpic_hydro() {
    name      = "vpic_hydro";
    nr_comp   = VPIC_HYDRO_N_COMP;
#if 0
    fld_names = { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
		  "px_nc", "py_nc", "pz_nc", "ke_nc",
		  "txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
		  "tzx_nc", "txy_nc", "_pad0", "_pad1", };
#endif
    run_all   = run_all_vpic_hydro;
    flags     = POFI_BY_KIND;
  }
} psc_output_fields_item_vpic_hydro_ops;

