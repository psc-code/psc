
#include "psc_output_fields_item_private.h"

#include <math.h>

#include "common_moments.cxx"

static void
run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	struct psc_mparticles *mprts_base, struct psc_mfields *mres_base,
	void (*do_run)(int p, fields_t flds, particle_range_t prts))
{
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mfields_t mf_res = mres_base->get_as<mfields_t>(0, 0);
  for (int p = 0; p < mprts.n_patches(); p++) {
    mf_res[p].zero();
    do_run(p, mf_res[p], mprts[p].range());
  }

  mf_res.put_as(mres_base, 0, mres_base->nr_fields);
  mprts.put_as(mprts_base, MP_DONT_COPY);
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
    int m = prt->kind();
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
    int m = prt->kind();
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
    int mm = prt->kind() * 3;

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



