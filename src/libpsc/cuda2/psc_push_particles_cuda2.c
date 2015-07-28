
#include "psc_push_particles_private.h"

#include "psc_cuda2.h"

#include <string.h>

// ======================================================================
// psc_push_particles: subclass "1vb_cuda2"

static void
psc_push_particles_1vbec_push_mprts_yz(struct psc_push_particles *push,
				     struct psc_mparticles *mprts_base,	
				     struct psc_mfields *mflds_base)					
{									
  // FIXME
  struct psc_mparticles *mprts =
    psc_mparticles_get_as(mprts_base, "cuda2", 0);
  struct psc_mfields *mflds =
    psc_mfields_get_as(mflds_base, "cuda2", EX, EX + 6);

  cuda2_1vbec_push_mprts_yz(mprts, mflds);

  psc_mparticles_put_as(mprts, mprts_base, 0);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);	
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vbec_cuda2"
									
struct psc_push_particles_ops						
psc_push_particles_1vbec_cuda2_ops = {	
  .name                  = "1vbec_cuda2",
  .push_mprts_yz         = psc_push_particles_1vbec_push_mprts_yz,
};
