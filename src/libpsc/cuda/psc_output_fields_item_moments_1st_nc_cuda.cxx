
#include "psc_output_fields_item_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "fields_item.hxx"

// ======================================================================
// rho_1st_nc_cuda

struct FieldsItem_rho_1st_nc_cuda : FieldsItemBase
{
  static const char* name() { return "rho_1st_nc_cuda"; }
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "rho_nc_cuda" }; } // FIXME
  constexpr static int flags = POFI_ADD_GHOSTS;

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base,
	   PscMfieldsBase mres_base) override
  {
    PscMparticlesCuda mprts = mprts_base.get_as<PscMparticlesCuda>();
    PscMfieldsCuda mf_res = mres_base.get_as<PscMfieldsCuda>(0, 0);
    cuda_mparticles *cmprts = mprts->cmprts();
    cuda_mfields *cmres = mf_res->cmflds;
    
    mf_res->zero();

    cuda_moments_yz_rho_1st_nc(cmprts, cmres);
    
    mprts.put_as(mprts_base, MP_DONT_COPY);
    mf_res.put_as(mres_base, 0, 1);
  }
};

FieldsItemOps<FieldsItem_rho_1st_nc_cuda> psc_output_fields_item_rho_1st_nc_cuda_ops;

// ======================================================================
// n_1st_cuda

struct FieldsItem_n_1st_cuda : FieldsItemBase
{
  static const char* name() { return "n_1st_cuda"; }
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n_1st_cuda" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base,
	   PscMfieldsBase mres_base) override
  {
    PscMparticlesCuda mprts = mprts_base.get_as<PscMparticlesCuda>();
    PscMfieldsCuda mres = mres_base.get_as<PscMfieldsCuda>(0, 0);
    cuda_mparticles *cmprts = mprts->cmprts();
    cuda_mfields *cmres = mres->cmflds;
    
    mres->zero();
    
    cuda_moments_yz_n_1st(cmprts, cmres);
    
    mprts.put_as(mprts_base, MP_DONT_COPY);
    mres.put_as(mres_base, 0, 1);
  }
};

FieldsItemOps<FieldsItem_n_1st_cuda> psc_output_fields_item_n_1st_cuda_ops;
