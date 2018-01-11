
#include "psc_output_fields_item_private.h"
#include "psc_fields_cuda.h"

#include "cuda_iface.h"

// ======================================================================

static void
calc_dive_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres_base)
{
  assert(ppsc->domain.gdims[0] == 1);

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX+3);
  mfields_cuda_t mf_res = mres_base->get_as<mfields_cuda_t>(0, 0);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mf.mflds())->cmflds;
  struct cuda_mfields *cmres = psc_mfields_cuda(mf_res.mflds())->cmflds;

  for (int p = 0; p < mf_res.n_patches(); p++) {
    cuda_mfields_calc_dive_yz(cmflds, cmres, p);
  }

  mf.put_as(mflds_base, 0, 0);
  mf_res.put_as(mres_base, 0, 1);
}

struct psc_output_fields_item_ops_dive_cuda : psc_output_fields_item_ops {
  psc_output_fields_item_ops_dive_cuda() {
    name         = "dive_cuda";
    nr_comp      = 1;
    fld_names[0] = "dive";
    run_all      = calc_dive_nc;
  }
} psc_output_fields_item_dive_cuda_ops;

