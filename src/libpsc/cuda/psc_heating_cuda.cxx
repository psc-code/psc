
#include <psc_heating_private.h>

#include "cuda_iface.h"
#include "psc_particles_cuda.h"

#include <stdlib.h>
#include <string.h>

// ======================================================================
// psc_heating subclass "cuda"

// ----------------------------------------------------------------------
// psc_heating_cuda_setup

static void
psc_heating_cuda_setup(struct psc_heating *heating)
{
  psc_heating_setup_super(heating);

  struct cuda_heating_foil foil;
  double val;

  psc_heating_spot_get_param_double(heating->spot, "zl", &val);
  foil.zl = val;
  psc_heating_spot_get_param_double(heating->spot, "zh", &val);
  foil.zh = val;
  psc_heating_spot_get_param_double(heating->spot, "xc", &val);
  foil.xc = val;
  psc_heating_spot_get_param_double(heating->spot, "yc", &val);
  foil.yc = val;
  psc_heating_spot_get_param_double(heating->spot, "rH", &val);
  foil.rH = val;
  psc_heating_spot_get_param_double(heating->spot, "T", &val);
  foil.T = val;
  psc_heating_spot_get_param_double(heating->spot, "Mi", &val);
  foil.Mi = val;
  foil.kind = heating->kind;

  foil.heating_dt = heating->every_step * ppsc->dt;

  cuda_heating_setup_foil(&foil);
}

// ----------------------------------------------------------------------
// psc_heating_cuda_run

void
psc_heating_cuda_run(struct psc_heating *heating, struct psc_mparticles *mprts_base,
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
  // FIXME, the above could be done by the superclass

  assert(strcmp(psc_mparticles_type(mprts_base), "cuda") == 0);
  struct psc_mparticles *mprts = mprts_base;

  struct cuda_mparticles *cmprts = mparticles_cuda_t(mprts).sub_;
  assert(strcmp(psc_heating_spot_type(heating->spot), "foil") == 0);
  cuda_heating_run_foil(cmprts);
}

// ----------------------------------------------------------------------
// psc_heating "cuda"

struct psc_heating_ops_cuda : psc_heating_ops {
  psc_heating_ops_cuda() {
    name                = "cuda";
    setup               = psc_heating_cuda_setup;
    run                 = psc_heating_cuda_run;
  }
} psc_heating_ops_cuda;
