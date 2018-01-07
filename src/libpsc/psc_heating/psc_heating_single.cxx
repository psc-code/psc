
#include <psc_heating_private.h>

#include <psc_particles_as_single.h>

#include <stdlib.h>

// ======================================================================
// psc_heating subclass "single"

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
// psc_heating_single_run

void
psc_heating_single_run(struct psc_heating *heating, struct psc_mparticles *mprts_base,
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
    particle_range_t prts = particle_range_mprts(mprts, p);
    struct psc_patch *patch = &psc->patch[p];
    PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
      particle_t *prt = particle_iter_deref(prt_iter);
      if (particle_kind(prt) != kind) {
	continue;
      }
      
      double xx[3] = {
	(&particle_x(prt))[0] + patch->xb[0],
	(&particle_x(prt))[1] + patch->xb[1],
	(&particle_x(prt))[2] + patch->xb[2],
      };

      double H = psc_heating_spot_get_H(heating->spot, xx);
      if (H > 0) {
	psc_heating_particle_kick(heating, prt, H);
      }
    }
  }

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ----------------------------------------------------------------------
// psc_heating "single"

struct psc_heating_ops psc_heating_ops_single = {
  .name                = "single",
  .run                 = psc_heating_single_run,
};

