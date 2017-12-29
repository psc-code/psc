
#include "psc_output_fields_item_private.h"

#include <math.h>

#include "common_moments.c"

static void
run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	struct psc_mparticles *mprts_base, struct psc_mfields *mres_base,
	void (*do_run)(int p, fields_t flds, particle_range_t prts))
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  struct psc_mfields *mres = psc_mfields_get_as(mres_base, FIELDS_TYPE, 0, 0);

  psc_mfields_view(mres);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t res = fields_t_mflds(mres, p);
    fields_t_zero_range(res, 0, mres->nr_fields);
    do_run(p, res, particle_range_mprts(mprts, p));
  }

  psc_mfields_put_as(mres, mres_base, 0, mres->nr_fields);
  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
}

// ======================================================================
// n

static void
do_n_run(int p, fields_t flds, particle_range_t prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int m = particle_kind(prt);
    DEPOSIT_TO_GRID_1ST_NC(prt, flds, m, 1.f);
  }
}

static void
n_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	  struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_n_run);
}

// ======================================================================
// rho

static void
do_rho_run(int p, fields_t flds, particle_range_t prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int m = particle_kind(prt);
    DEPOSIT_TO_GRID_1ST_NC(prt, flds, 0, ppsc->kinds[m].q);
  }
}

static void
rho_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	  struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_rho_run);
}

// ======================================================================
// v

static void
do_v_run(int p, fields_t flds, particle_range_t prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 3;

    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_NC(prt, flds, mm + m, vxi[m]);
    }
  }
}

static void
v_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	  struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_v_run);
}



