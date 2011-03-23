
#include "psc.h"

#include <stdlib.h>

// ----------------------------------------------------------------------
// mfields_base_alloc

void
mfields_base_alloc(mfields_base_t *flds, int nr_fields)
{
  flds->f = calloc(psc.nr_patches, sizeof(*flds->f));
  foreach_patch(p) {
    int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
    int ihg[3] = { psc.patch[p].ldims[0] + psc.ibn[0],
		   psc.patch[p].ldims[1] + psc.ibn[1],
		   psc.patch[p].ldims[2] + psc.ibn[2] };
    fields_base_alloc(&flds->f[p], ilg, ihg, nr_fields);
  }
}

// ----------------------------------------------------------------------
// mfields_base_destroy

void
mfields_base_destroy(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_free(&flds->f[p]);
  }
  free(flds->f);
}

