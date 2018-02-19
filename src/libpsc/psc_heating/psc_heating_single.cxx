
#include <psc_heating_private.h>

#include <psc_particles_as_single.h>

using real_t = mparticles_t::real_t;

#include <stdlib.h>

// ======================================================================
// psc_heating subclass "single"

// ----------------------------------------------------------------------
// psc_heating_particle_kick

static void
psc_heating_particle_kick(struct psc_heating *heating, particle_t *prt, real_t H)
{
  struct psc *psc = ppsc;

  real_t heating_dt = heating->every_step * psc->dt;

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

  real_t ranx = sqrtf(-2.f*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
  real_t rany = sqrtf(-2.f*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4);
  real_t ranz = sqrtf(-2.f*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6);

  real_t Dpxi = sqrtf(H * heating_dt);
  real_t Dpyi = sqrtf(H * heating_dt);
  real_t Dpzi = sqrtf(H * heating_dt);

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

  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  
  psc_foreach_patch(psc, p) {
    mparticles_t::patch_t& prts = mprts[p];
    Grid_t::Patch& patch = psc->grid().patches[p];
    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      particle_t *prt = &*prt_iter;
      if (prt->kind() != kind) {
	continue;
      }
      
      double xx[3] = {
	prt->xi + patch.xb[0],
	prt->yi + patch.xb[1],
	prt->zi + patch.xb[2],
      };

      double H = psc_heating_spot_get_H(heating->spot, xx);
      if (H > 0) {
	psc_heating_particle_kick(heating, prt, H);
      }
    }
  }

  mprts.put_as(mprts_base);
}

// ----------------------------------------------------------------------
// psc_heating "single"

struct psc_heating_ops_single : psc_heating_ops {
  psc_heating_ops_single() {
    name                = "single";
    run                 = psc_heating_single_run;
  }
} psc_heating_ops_single;

