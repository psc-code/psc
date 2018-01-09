
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_event_generator, struct psc_event_generator);

void psc_event_generator_run(struct psc_event_generator *gen,
			     struct psc_mparticles *mprts, struct psc_mfields *mflds);
