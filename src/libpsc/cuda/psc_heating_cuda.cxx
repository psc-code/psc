
#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "heating.hxx"

#include <stdlib.h>
#include <string.h>

// ======================================================================
// psc_heating subclass "cuda"

struct HeatingCuda : HeatingBase
{
  // ----------------------------------------------------------------------
  // ctor

  template<typename FUNC>
  HeatingCuda(int every_step, int tb, int te, int kind, FUNC get_H)
    : every_step_(every_step),
      tb_(tb), te_(te),
      kind_(kind)
  {
    struct cuda_heating_foil foil;
    double val;

#if 0
    assert(strcmp(psc_heating_spot_type(&spot_), "foil") == 0);
    psc_heating_spot_get_param_double(&spot_, "zl", &val);
    foil.zl = val;
    psc_heating_spot_get_param_double(&spot_, "zh", &val);
    foil.zh = val;
    psc_heating_spot_get_param_double(&spot_, "xc", &val);
    foil.xc = val;
    psc_heating_spot_get_param_double(&spot_, "yc", &val);
    foil.yc = val;
    psc_heating_spot_get_param_double(&spot_, "rH", &val);
    foil.rH = val;
    psc_heating_spot_get_param_double(&spot_, "T", &val);
    foil.T = val;
    psc_heating_spot_get_param_double(&spot_, "Mi", &val);
    foil.Mi = val;
    foil.kind = kind_;
    foil.heating_dt = every_step_ * ppsc->dt;
#endif
    cuda_heating_setup_foil(&foil);
  }

  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base) override
  {
    struct psc *psc = ppsc;

    // only heating between heating_tb and heating_te
    if (psc->timestep < tb_ || psc->timestep >= te_) {
      return;
    }

    if (psc->timestep % every_step_ != 0) {
      return;
    }

    PscMparticlesCuda mprts = mprts_base.get_as<PscMparticlesCuda>();
    cuda_mparticles *cmprts = mprts->cmprts();
    cuda_heating_run_foil(cmprts);
  }

private:
  int every_step_;
  int tb_, te_;
  int kind_;
};

