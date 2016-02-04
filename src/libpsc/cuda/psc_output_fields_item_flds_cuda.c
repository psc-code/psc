
#include "psc_output_fields_item_private.h"

#include "psc_cuda.h"

// ======================================================================

static void
calc_dive_nc(struct psc_output_fields_item *item, struct psc_fields *flds_base,
	     struct psc_particles *prts, struct psc_fields *f_base)
{
  assert(ppsc->domain.gdims[0] == 1);

  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, EX + 3);
  struct psc_fields *f = psc_fields_get_as(f_base, "cuda", 0, 0);

  cuda_calc_dive_yz(flds, f);

  psc_fields_put_as(flds, flds_base, 0, 0);
  psc_fields_put_as(f, f_base, 0, 1);
}

struct psc_output_fields_item_ops psc_output_fields_item_dive_cuda_ops = {
  .name      = "dive_cuda",
  .nr_comp   = 1,
  .fld_names = { "dive" },
  .run       = calc_dive_nc,
};

