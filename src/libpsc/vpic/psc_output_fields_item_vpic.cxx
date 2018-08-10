
#include "psc_output_fields_item_private.h"

#include "vpic_config.h"
#include "fields.hxx"
#include "fields_item.hxx"
#include "psc_fields_as_single.h"
#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"

#include "fields_item_vpic.hxx"

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

  Moment_vpic_hydro(const Grid_t& grid, MPI_Comm comm)
    : Base(grid, comm)
  {}
  
  void run(MparticlesVpic& mprts)
  {
    const auto& kinds = mprts.grid().kinds;
    auto& mres = this->mres_;
    auto& grid = mprts.grid();
    auto mf_hydro = MfieldsVpic{ppsc->grid(), 16, { 1, 1, 1 }};
    Simulation *sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    for (int kind = 0; kind < kinds.size(); kind++) {
      HydroArray *vmflds_hydro = mf_hydro.vmflds_hydro;
      sim->moments_run(vmflds_hydro, &mprts.vmprts_, kind);
      
      auto& mf = mf_hydro.get_as<Mfields>(0, VPIC_HYDRO_N_COMP);
      for (int p = 0; p < mres.n_patches(); p++) {
	Fields F(mf[p]), R(mres[p]);
	grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
	    for (int m = 0; m < VPIC_HYDRO_N_COMP; m++) {
	      R(m + kind * VPIC_HYDRO_N_COMP, ix,iy,iz) = F(m, ix,iy,iz);
	    }
	  });
      }
      mf_hydro.put_as(mf, 0, 0);
    }
  }
};

FieldsItemOps<FieldsItemMoment<Moment_vpic_hydro>> psc_output_fields_item_vpic_hydro_ops;

