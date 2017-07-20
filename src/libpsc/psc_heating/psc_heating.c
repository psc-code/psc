
#include "psc_heating_private.h"

#include <psc_particles_as_single.h> // FIXME

#include <stdlib.h>

// ======================================================================
// psc_heating

// ----------------------------------------------------------------------
// _psc_heating_setup

static void
_psc_heating_setup(struct psc_heating *heating)
{
  double width = heating->zh - heating->zl;
  heating->fac = (8.f * pow(heating->T, 1.5)) / (sqrt(heating->Mi) * width);
  // FIXME, I don't understand the sqrt(Mi) in here
}

// ----------------------------------------------------------------------
// psc_heating_get_H

static particle_real_t
psc_heating_get_H(struct psc_heating *heating, particle_real_t *xx)
{
  particle_real_t zl = heating->zl;
  particle_real_t zh = heating->zh;
  particle_real_t xc = heating->xc;
  particle_real_t yc = heating->yc;
  particle_real_t rH = heating->rH;
  particle_real_t fac = heating->fac;
  particle_real_t x = xx[0], y = xx[1], z = xx[2];

  if (z <= zl || z >= zh) {
    return 0;
  }

  return fac * exp(-(sqr(x-xc) + sqr(y-yc)) / sqr(rH));
}
  
// ----------------------------------------------------------------------
// psc_heating_particle_kick

static void
psc_heating_particle_kick(struct psc_heating *heating, particle_t *prt, particle_real_t H)
{
  struct psc *psc = ppsc;

  particle_real_t heating_dt = heating->every_step * psc->dt;

  float ran1, ran2, ran3, ran4, ran5, ran6;
  do {
    ran1 = random() / ((float) RAND_MAX + 1);
    ran2 = random() / ((float) RAND_MAX + 1);
    ran3 = random() / ((float) RAND_MAX + 1);
    ran4 = random() / ((float) RAND_MAX + 1);
    ran5 = random() / ((float) RAND_MAX + 1);
    ran6 = random() / ((float) RAND_MAX + 1);
  } while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
	   ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);

  particle_real_t ranx = sqrtf(-2.f*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
  particle_real_t rany = sqrtf(-2.f*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4);
  particle_real_t ranz = sqrtf(-2.f*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6);

  particle_real_t Dpxi = sqrtf(H * heating_dt);
  particle_real_t Dpyi = sqrtf(H * heating_dt);
  particle_real_t Dpzi = sqrtf(H * heating_dt);

  prt->pxi += Dpxi * ranx;
  prt->pyi += Dpyi * rany;
  prt->pzi += Dpzi * ranz;
}

// ----------------------------------------------------------------------
// psc_heating_run

void
psc_heating_run(struct psc_heating *heating, struct psc_mparticles *mprts_base,
		struct psc_mfields *mflds_base)
{
  struct psc *psc = ppsc;

  // only heating between heating_tb and heating_te
  if (psc->timestep < heating->tb || psc->timestep >= heating->te) {
    return;
  }

  if (psc->timestep % heating->every_step != 0) {
    return;
  }

  int kind = heating->kind;

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  
  psc_foreach_patch(psc, p) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_patch *patch = &psc->patch[p];
    for (int n = 0; n < prts->n_part; n++) {
      particle_t *prt = particles_get_one(prts, n++);
      if (particle_kind(prt) != kind) {
	continue;
      }
      
      particle_real_t xx[3] = {
	(&particle_x(prt))[0] + patch->xb[0],
	(&particle_x(prt))[1] + patch->xb[1],
	(&particle_x(prt))[2] + patch->xb[2],
      };

      double H = psc_heating_get_H(heating, xx);
      if (H > 0) {
	psc_heating_particle_kick(heating, prt, H);
      }
    }
  }

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ----------------------------------------------------------------------
// psc_heating class

struct mrc_class_psc_heating mrc_class_psc_heating = {
  .name             = "psc_heating",
  .size             = sizeof(struct psc_heating),
  .param_descr      = psc_heating_descr,
  .setup            = _psc_heating_setup,
};

