
#include "psc_method_private.h"

#include <psc_fields_vpic.h>
#include <psc_particles_vpic.h>
#include <psc_push_particles_vpic.h>

#include <vpic_iface.h>


// ======================================================================
// psc_method "vpic"

static void
psc_method_vpic_initialize(struct psc_method *method, struct psc *psc)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(psc->flds)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(psc->particles)->vmprts;
  struct vpic_push_particles *vpushp = psc_push_particles_vpic(psc->push_particles)->vpushp;

  vpic_simulation_init2(vpushp, vmflds, vmprts);

  mpi_printf(psc_comm(psc), "Initialization complete.\n");
}

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops psc_method_ops_vpic = {
  .name                = "vpic",
  .initialize          = psc_method_vpic_initialize,
};
