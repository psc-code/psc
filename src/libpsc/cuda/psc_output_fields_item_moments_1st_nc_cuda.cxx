
#include "psc_output_fields_item_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

// ======================================================================
// rho_1st_nc

static void
rho_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	    struct psc_mparticles *mprts_base, struct psc_mfields *mres_base)
{
  mparticles_cuda_t mprts = mprts_base->get_as<mparticles_cuda_t>();
  mfields_cuda_t mf_res = mres_base->get_as<mfields_cuda_t>(0, 0);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts.mprts())->cmprts;
  struct cuda_mfields *cmres = psc_mfields_cuda(mf_res.mflds())->cmflds;
    
  mf_res.zero();

  cuda_moments_yz_rho_1st_nc(cmprts, cmres);

  mprts.put_as(mprts_base, MP_DONT_COPY);
  mf_res.put_as(mres_base, 0, 1);
}

// ----------------------------------------------------------------------
// psc_output_fields_item: subclass "rho_1st_nc_cuda"

struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_cuda_ops = {
  .name               = "rho_1st_nc_cuda",
  .nr_comp            = 1,
  .fld_names          = { "rho_nc_cuda" }, // FIXME
  .run_all            = rho_run_all,
  .flags              = POFI_ADD_GHOSTS,
};

// ======================================================================
// n_1st

// ----------------------------------------------------------------------
// n_1st_run_all

static void
n_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts_base, struct psc_mfields *mres_base)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "cuda", 0);
  mfields_cuda_t mf_res = mres_base->get_as<mfields_cuda_t>(0, 0);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  struct cuda_mfields *cmres = psc_mfields_cuda(mf_res.mflds())->cmflds;

  mf_res.zero();

  cuda_moments_yz_n_1st(cmprts, cmres);

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
  mf_res.put_as(mres_base, 0, 1);
}

// ----------------------------------------------------------------------
// psc_output_fields_item: subclass "n_1st_cuda"

struct psc_output_fields_item_ops psc_output_fields_item_n_1st_cuda_ops = {
  .name               = "n_1st_cuda",
  .nr_comp            = 1,
  .fld_names          = { "n" },
  .run_all            = n_1st_run_all,
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,
};

