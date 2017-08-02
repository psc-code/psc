
#include "psc_output_fields_item_private.h"

#include "psc_cuda.h"

// ======================================================================

static void
calc_dive_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres_base)
{
  assert(ppsc->domain.gdims[0] == 1);

  for (int p = 0; p < mres_base->nr_patches; p++) {
    struct psc_fields *flds_base = psc_mfields_get_patch(mflds_base, p);
    struct psc_fields *res_base = psc_mfields_get_patch(mres_base, p);
    struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, EX + 3);
    struct psc_fields *res = psc_fields_get_as(res_base, "cuda", 0, 0);

    cuda_calc_dive_yz(flds, res);

    psc_fields_put_as(flds, flds_base, 0, 0);
    psc_fields_put_as(res, res_base, 0, 1);
  }
}

struct psc_output_fields_item_ops psc_output_fields_item_dive_cuda_ops = {
  .name      = "dive_cuda",
  .nr_comp   = 1,
  .fld_names = { "dive" },
  .run_all   = calc_dive_nc,
};

