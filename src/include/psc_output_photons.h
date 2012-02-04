
#include <mrc_obj.h>

#include "psc.h"
#include "psc_photons.h"

MRC_CLASS_DECLARE(psc_output_photons, struct psc_output_photons);

void psc_output_photons_run(struct psc_output_photons *output_photons,
			      mphotons_t *photons);
