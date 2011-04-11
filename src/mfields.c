
#include "psc.h"

#include <stdlib.h>

// ----------------------------------------------------------------------
// mfields_base_alloc

mfields_base_t *
mfields_base_alloc(struct mrc_domain *domain,
		   int nr_fields, int ibn[3])
{
  mfields_base_t *flds = calloc(1, sizeof(*flds));
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &flds->nr_patches);
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));
  for (int p = 0; p < flds->nr_patches; p++) {
    int ilg[3] = { -ibn[0], -ibn[1], -ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + ibn[0],
		   patches[p].ldims[1] + ibn[1],
		   patches[p].ldims[2] + ibn[2] };
    fields_base_alloc(&flds->f[p], ilg, ihg, nr_fields);
  }
  return flds;
}

// ----------------------------------------------------------------------
// mfields_base_destroy

void
mfields_base_destroy(mfields_base_t *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_base_free(&flds->f[p]);
  }
  free(flds->f);
  free(flds);
}

