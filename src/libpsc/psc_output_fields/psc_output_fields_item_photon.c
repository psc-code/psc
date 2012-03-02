
#include "psc_output_fields_item_private.h"

#include "psc_moments.h"

// ======================================================================

static void
calc_photon_n(struct psc_output_fields_item *item, mfields_base_t *flds,
	      mparticles_base_t *particles, mfields_c_t *res)
{
  psc_moments_calc_photon_n(ppsc->moments, ppsc->mphotons, res);
}

struct psc_output_fields_item_ops psc_output_fields_item_photon_n_ops = {
  .name      = "photon_n",
  .nr_comp   = 1,
  .fld_names = { "photon_n" },
  .run       = calc_photon_n,
};

