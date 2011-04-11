
#include "psc.h"

#include <math.h>

static void
push_photons_patch(int p, photons_t *photons)
{
  for (int i = 0; i < photons->nr; i++) {
    photon_t *p = photons_get_one(photons, i);

    real p_tot = sqrt(sqr(p->p[0]) + sqr(p->p[1]) + sqr(p->p[2]));
    
    for (int d = 0; d < 3; d++) {
      real v = p->p[d] / p_tot;
      p->x[d] += v * psc.dt;
    }
  }
}

void
psc_push_photons_run(mphotons_t *mphotons)
{
  psc_foreach_patch(&psc, p) {
    push_photons_patch(p, &mphotons->p[p]);
  }
}

