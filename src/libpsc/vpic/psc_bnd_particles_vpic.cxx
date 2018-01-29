
#include "psc_bnd_particles_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_bnd_particles_vpic_exchange_particles

static void
psc_bnd_particles_vpic_exchange_particles(struct psc_bnd_particles *bnd,
					  struct psc_mparticles *mprts_base)
{
  // nothing to be done here since the particle exchange happens in
  // push_particles_run() already, at least for now...
}

// ----------------------------------------------------------------------
// psc_bnd_particles_vpic_reset

static void
psc_bnd_particles_vpic_reset(struct psc_bnd_particles *bnd)
{
}

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "vpic"

struct psc_bnd_particles_ops_vpic : psc_bnd_particles_ops {
  psc_bnd_particles_ops_vpic() {
    name                  = "vpic";
    exchange_particles    = psc_bnd_particles_vpic_exchange_particles;
    reset                 = psc_bnd_particles_vpic_reset;
  }
} psc_bnd_particles_vpic_ops;

