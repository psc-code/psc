
#include "psc_output_fields_item_private.h"

#include "psc_cuda.h"

// ======================================================================

static void
calc_dive_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres_base)
{
  assert(ppsc->domain.gdims[0] == 1);

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", EX, EX + 3);
  struct psc_mfields *mres = psc_mfields_get_as(mres_base, "cuda", 0, 0);
  for (int p = 0; p < mres->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields *res = psc_mfields_get_patch(mres, p);

    cuda_calc_dive_yz(flds, res);
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
  psc_mfields_put_as(mres, mres_base, 0, 1);
}

struct psc_output_fields_item_ops psc_output_fields_item_dive_cuda_ops = {
  .name      = "dive_cuda",
  .nr_comp   = 1,
  .fld_names = { "dive" },
  .run_all   = calc_dive_nc,
};

