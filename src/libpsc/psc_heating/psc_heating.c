
#include "psc_heating_private.h"

#include <psc_particles_as_single.h> // FIXME

#include <mrc_profile.h>
#include <stdlib.h>

// ======================================================================
// psc_heating

// ----------------------------------------------------------------------
// psc_heating_get_spot

struct psc_heating_spot *
psc_heating_get_spot(struct psc_heating *heating)
{
  return heating->spot;
}

// ----------------------------------------------------------------------
// psc_heating_run

void
psc_heating_run(struct psc_heating *heating, struct psc_mparticles *mprts_base,
		struct psc_mfields *mflds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("heating_run", 1., 0, 0);
  }  

  prof_start(pr);
  struct psc_heating_ops *ops = psc_heating_ops(heating);

  assert(ops && ops->run);
  ops->run(heating, mprts_base, mflds_base);

  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_heating_init

static void
psc_heating_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_heating, &psc_heating_ops_single);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_heating, &psc_heating_ops_cuda);
#endif
}

// ----------------------------------------------------------------------
// psc_heating class

struct mrc_class_psc_heating mrc_class_psc_heating = {
  .name             = "psc_heating",
  .size             = sizeof(struct psc_heating),
  .param_descr      = psc_heating_descr,
  .init             = psc_heating_init,
};



