
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
// psc_bnd_particles_vpic_unsetup

static void
psc_bnd_particles_vpic_unsetup(struct psc_bnd_particles *bnd)
{
}

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "vpic"

struct psc_bnd_particles_ops psc_bnd_particles_vpic_ops = {
  .name                  = "vpic",
  .exchange_particles    = psc_bnd_particles_vpic_exchange_particles,
  .unsetup               = psc_bnd_particles_vpic_unsetup,
};

