
#include "psc_event_generator_private.h"

#include <mrc_profile.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// event_generator_run_patch
//
// create new photons on given patch

static void
event_generator_run_patch(int p, photons_t *photons)
{
  for (;;) { // loop forever
    // until the random number is > .99, that is, so it'll loop about 100 times
    // on average
    float r = random() / (float) RAND_MAX;
    if (r > .99)
      break;

    photons_realloc(photons, photons->nr + 1);
    photon_t *ph = photons_get_one(photons, photons->nr++);
    ph->x[0] = 5. * 1e-6 / ppsc->coeff.ld;
    ph->x[1] = 5. * 1e-6 / ppsc->coeff.ld;
    ph->x[2] = 5. * 1e-6 / ppsc->coeff.ld;
    ph->p[0] = 0.1;
    ph->p[1] = 0.;
    ph->p[2] = 0.;
    ph->wni = 1.;
  }
}

// ----------------------------------------------------------------------
// psc_event_generator_demo_run

static void
psc_event_generator_demo_run(struct psc_event_generator *gen,
			     struct psc_mparticles *mparticles_base,
			  mfields_base_t *mflds_base,
			  mphotons_t *mphotons)
{
  psc_foreach_patch(ppsc, p) {
    event_generator_run_patch(p, &mphotons->p[p]);
  }
}

// ======================================================================
// psc_event_generator: subclass "demo"

struct psc_event_generator_ops psc_event_generator_demo_ops = {
  .name                  = "demo",
  .run                   = psc_event_generator_demo_run,
};
