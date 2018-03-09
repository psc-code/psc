
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
// psc_heating_init

extern struct psc_heating_ops psc_heating_ops_single;
extern struct psc_heating_ops psc_heating_ops_double;
extern struct psc_heating_ops psc_heating_ops_cuda;

static void
psc_heating_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_heating, &psc_heating_ops_single);
  mrc_class_register_subclass(&mrc_class_psc_heating, &psc_heating_ops_double);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_heating, &psc_heating_ops_cuda);
#endif
}

// ----------------------------------------------------------------------
// psc_heating class

struct mrc_class_psc_heating_ : mrc_class_psc_heating {
  mrc_class_psc_heating_() {
    name             = "psc_heating";
    size             = sizeof(struct psc_heating);
    param_descr      = psc_heating_descr;
    init             = psc_heating_init;
  }
} mrc_class_psc_heating;



