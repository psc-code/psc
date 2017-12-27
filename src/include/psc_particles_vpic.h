
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct psc_mparticles_vpic {
  Particles *vmprts;
  Simulation *sim;
};

#define psc_mparticles_vpic(mprts)({					\
      assert((struct psc_mparticles_ops *) mprts->obj.ops == &psc_mparticles_vpic_ops); \
      mrc_to_subobj(mprts, struct psc_mparticles_vpic);			\
})

#endif
