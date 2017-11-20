
#include "psc_output_fields_item_private.h"

#include "psc_fields_as_c.h"
#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"
#include "vpic_iface.h"


// ----------------------------------------------------------------------
// run_all_vpic_fields

static void
run_all_vpic_fields(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		   struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, 0, 16);
  
  for (int p = 0; p < mres->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    fields_t res = fields_t_mflds(mres, p);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      for (int m = 0; m < 16; m++) {
	_F3(res, m, ix,iy,iz) = _F3(flds, m, ix,iy,iz);
      }
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_output_fields_item "vpic_fields"

struct psc_output_fields_item_ops psc_output_fields_item_vpic_fields_ops = {
  .name      = "vpic_fields",
  .nr_comp   = 16,
  .fld_names = { "ex_ec", "ey_ec", "ez_ec", "dive_nc",
                 "hx_fc", "hy_fc", "hz_fc", "divb_cc",
                 "tcax", "tcay", "tcaz", "rhob_nc",
                 "jx_ec", "jy_ec", "jz_ec", "rho_nc", },
  .run_all   = run_all_vpic_fields,
};

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
  psc_mfields_set_param_int3(mflds_hydro, "ibn", (int [3]) { 1, 1, 1});
  psc_mfields_setup(mflds_hydro);

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "vpic", 0);

  Simulation *sim;
  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);

  for (int kind = 0; kind < ppsc->nr_kinds; kind++) {
    struct HydroArray *vmflds_hydro = psc_mfields_vpic(mflds_hydro)->vmflds_hydro;
    Particles *vmprts = psc_mparticles_vpic(mprts)->vmprts;
    Simulation_moments_run(sim, vmflds_hydro, vmprts, kind);
    
    struct psc_mfields *mflds = psc_mfields_get_as(mflds_hydro, FIELDS_TYPE, 0, VPIC_HYDRO_N_COMP);
  
    for (int p = 0; p < mres->nr_patches; p++) {
      fields_t flds = fields_t_mflds(mflds, p);
      fields_t res = fields_t_mflds(mres, p);
      psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	for (int m = 0; m < VPIC_HYDRO_N_COMP; m++) {
	  _F3(res, m + kind * VPIC_HYDRO_N_COMP, ix,iy,iz) = _F3(flds, m, ix,iy,iz);
	}
      } foreach_3d_end;
    }
    psc_mfields_put_as(mflds, mflds_hydro, 0, 0);
  }

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
  
  psc_mfields_destroy(mflds_hydro);
}

// ----------------------------------------------------------------------
// psc_output_fields_item "vpic_hydro"

struct psc_output_fields_item_ops psc_output_fields_item_vpic_hydro_ops = {
  .name      = "vpic_hydro",
  .nr_comp   = VPIC_HYDRO_N_COMP,
  .fld_names = { "jx_nc", "jy_nc", "jz_nc", "rho_nc",
                 "px_nc", "py_nc", "pz_nc", "ke_nc",
                 "txx_nc", "tyy_nc", "tzz_nc", "tyz_nc",
                 "tzx_nc", "txy_nc", "_pad0", "_pad1", },
  .run_all   = run_all_vpic_hydro,
  .flags     = POFI_BY_KIND,
};

