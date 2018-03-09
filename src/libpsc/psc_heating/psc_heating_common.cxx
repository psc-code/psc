
#include <psc_heating_private.h>
#include "heating.hxx"

using real_t = mparticles_t::real_t;

#include <stdlib.h>

// ======================================================================
// psc_heating subclass "sub"

template<typename MP>
struct Heating_ : HeatingBase
{
  using Self = Heating_<MP>;
  
  // ----------------------------------------------------------------------
  // ctor

  Heating_(int every_step, int tb, int te, int kind,
	   psc_heating_spot& spot)
    : every_step_(every_step),
      tb_(tb), te_(te),
      kind_(kind),
      spot_(spot)
  {}
  
  // ----------------------------------------------------------------------
  // kick_particle

  void kick_particle(particle_t *prt, real_t H)
  {
    struct psc *psc = ppsc;

    real_t heating_dt = every_step_ * psc->dt;

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
  // run

  static void run(struct psc_heating *_heating, struct psc_mparticles *mprts_base,
		  struct psc_mfields *mflds_base)
  {
    PscHeating<Self> heating(_heating);
    struct psc *psc = ppsc;

    // only heating between heating_tb and heating_te
    if (psc->timestep < heating->tb_ || psc->timestep >= heating->te_) {
      return;
    }

    if (psc->timestep % heating->every_step_ != 0) {
      return;
    }

    int kind = heating->kind_;

    mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  
    psc_foreach_patch(psc, p) {
      mparticles_t::patch_t& prts = mprts[p];
      auto& patch = psc->grid().patches[p];
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

	double H = psc_heating_spot_get_H(&heating->spot_, xx);
	if (H > 0) {
	  heating->kick_particle(prt, H);
	}
      }
    }

    mprts.put_as(mprts_base);
  }

private:
  int every_step_;
  int tb_, te_;
  int kind_;
  psc_heating_spot& spot_;
};

// ----------------------------------------------------------------------
// psc_heating "sub"

struct psc_heating_ops_sub : psc_heating_ops {
  using Heating_t = Heating_<mparticles_t::sub_t>;
  using PscHeating_t = PscHeatingWrapper<Heating_t>;
  psc_heating_ops_sub() {
    name                = PARTICLE_TYPE;
    size                = PscHeating_t::size;
    setup               = PscHeating_t::setup;
    destroy             = PscHeating_t::destroy;
    run                 = Heating_t::run;
  }
} psc_heating_ops_sub;


