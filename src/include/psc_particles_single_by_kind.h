
#ifndef PSC_PARTICLES_SINGLE_BY_KIND_H
#define PSC_PARTICLES_SINGLE_BY_KIND_H

#include "psc_particles_private.h"

struct psc_mparticles_single_by_kind {
  bk_mparticles *bkmprts;
};

#define psc_mparticles_single_by_kind(mprts)({				\
      assert((struct psc_mparticles_ops *) mprts->obj.ops == &psc_mparticles_single_by_kind_ops); \
      mrc_to_subobj(mprts, struct psc_mparticles_single_by_kind);	\
})

#endif
