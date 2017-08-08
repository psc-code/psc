
#include "psc_output_photons_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_output_photons_none_run

static void
psc_output_photons_none_run(struct psc_output_photons *out,
				 mphotons_t *photons)
{
}

// ======================================================================
// psc_output_photons: subclass "none"

struct psc_output_photons_ops psc_output_photons_none_ops = {
  .name                  = "none",
  .run                   = psc_output_photons_none_run,
};
