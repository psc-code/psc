
#include "psc_output_fields_item_private.h"
#include "psc_bnd.h"

#include <math.h>

#include "common_moments.c"

// ======================================================================
// n

static void
do_n_run(int p, fields_t *pf, struct psc_particles *prts)
{
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / ppsc->dx[0], dyi = 1.f / ppsc->dx[1], dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_1ST_NC(part, pf, m, 1.f);
  }
}

static void
n_run(struct psc_output_fields_item *item, mfields_base_t *flds,
      mparticles_base_t *particles_base, mfields_c_t *res)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  psc_mfields_zero_range(res, 0, res->nr_fields);
  
  psc_foreach_patch(ppsc, p) {
    do_n_run(p, psc_mfields_get_patch(res, p),
	     psc_mparticles_get_patch(particles, p));
  }

  psc_mparticles_put_cf(particles, particles_base, MP_DONT_COPY);

  psc_bnd_add_ghosts(item->bnd, res, 0, res->nr_fields);
}

static int
n_get_nr_components(struct psc_output_fields_item *item)
{
  return ppsc->nr_kinds;
}

static const char *
n_get_component_name(struct psc_output_fields_item *item, int m)
{
  static char s[100];
  sprintf(s, "n_%s", ppsc->kinds[m].name);
  return s;
}

