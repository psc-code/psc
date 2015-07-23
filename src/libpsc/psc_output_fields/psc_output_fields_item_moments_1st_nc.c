
#include "psc_output_fields_item_private.h"

#include <math.h>

#include "common_moments.c"

// ======================================================================
// n

static void
do_n_run(int p, fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_1ST_NC(part, pf, m, 1.f);
  }
}

static void
n_run(struct psc_output_fields_item *item, struct psc_fields *flds,
      struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_n_run(res->p, res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

// ======================================================================
// rho

static void
do_rho_run(int p, fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_1ST_NC(part, pf, 0, ppsc->kinds[m].q);
  }
}

static void
rho_run_patches(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres_base)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  struct psc_mfields *mres = psc_mfields_get_as(mres_base, FIELDS_TYPE, 0, 0);
  psc_mfields_zero_range(mres, 0, mres->nr_fields);
  for (int p = 0; p < mres->nr_patches; p++) {
    do_rho_run(p, psc_mfields_get_patch(mres, p), psc_mparticles_get_patch(mprts, p));
  }
  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
  psc_mfields_put_as(mres, mres_base, 0, 1);
}

