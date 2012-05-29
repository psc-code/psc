
#include "psc.h"
#include "psc_particles_private.h"

#include <stdlib.h>
#include <assert.h>

// ----------------------------------------------------------------------
// psc_mfields_mix_setup

static void
psc_mfields_mix_setup(struct psc_mfields *mflds)
{
  mflds->comp_name = calloc(mflds->nr_fields, sizeof(*mflds->comp_name));

  struct mrc_patch *patches = mrc_domain_get_patches(mflds->domain,
						     &mflds->nr_patches);
  mflds->flds = calloc(mflds->nr_patches, sizeof(*mflds->flds));
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_fields_create(psc_mfields_comm(mflds));
    if ((p & 1) == 0) {
      psc_fields_set_type(flds, "single");
    } else {
      psc_fields_set_type(flds, "cuda");
    }
    for (int d = 0; d < 3; d++) {
      flds->ib[d] = -mflds->ibn[d];
      flds->im[d] = patches[p].ldims[d] + 2 * mflds->ibn[d];
    }
    flds->nr_comp = mflds->nr_fields;
    flds->first_comp = mflds->first_comp;
    flds->p = p;
    psc_fields_setup(flds);
    mflds->flds[p] = flds;
  }
}

// ======================================================================
// psc_mfields: subclass "mix"
  
struct psc_mfields_ops psc_mfields_mix_ops = {
  .name                    = "mix",
  .setup                   = psc_mfields_mix_setup,
};

