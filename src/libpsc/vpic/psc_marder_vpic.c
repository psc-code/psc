
#include "psc_marder_private.h"

#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_marder_vpic_run

static void
psc_marder_vpic_run(struct psc_marder *marder,
		    struct psc_mfields *mflds_base,
		    struct psc_mparticles *mprts_base)
{
  struct psc *psc = ppsc; // FIXME
  int step = psc->timestep;

  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts_base)->vmprts;
    
  // Divergence clean e
  int clean_div_e_interval = marder->clean_div_e_interval;
  if (clean_div_e_interval > 0 &&
      step % clean_div_e_interval == 0) {
    vpic_clean_div_e(vmflds, vmprts);
  }

  // Divergence clean b
  int clean_div_b_interval = marder->clean_div_b_interval;
  if (clean_div_b_interval > 0 &&
      step % clean_div_b_interval == 0) {
    vpic_clean_div_b(vmflds);
  }

  // Synchronize the shared faces
  int sync_shared_interval = marder->sync_shared_interval;
  if (sync_shared_interval > 0 &&
      step % sync_shared_interval == 0) {
    vpic_sync_faces(vmflds);
  }
}

// ----------------------------------------------------------------------
// psc_marder: subclass "vpic"

struct psc_marder_ops psc_marder_vpic_ops = {
  .name                  = "vpic",
  .run                   = psc_marder_vpic_run,
};

