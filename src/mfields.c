
#include "psc.h"

#include <stdlib.h>

// ----------------------------------------------------------------------
// mfields_base_alloc

void
mfields_base_alloc(struct mrc_domain *domain, mfields_base_t *flds,
		   int nr_fields, int ibn[3])
{
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &flds->nr_patches);
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));
  for (int p = 0; p < flds->nr_patches; p++) {
    int ilg[3] = { -ibn[0], -ibn[1], -ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + ibn[0],
		   patches[p].ldims[1] + ibn[1],
		   patches[p].ldims[2] + ibn[2] };
    fields_base_alloc(&flds->f[p], ilg, ihg, nr_fields);
  }
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
  flds->nr_patches = -1;
}

