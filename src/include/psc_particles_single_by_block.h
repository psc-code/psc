
#ifndef PSC_PARTICLE_SINGLE_BY_BLOCK_H
#define PSC_PARTICLE_SINGLE_BY_BLOCK_H

#include "psc_particles_private.h"

#include "particles.hxx"

#define PTYPE PTYPE_SINGLE_BY_BLOCK
#include "psc_particle_common.h"
using psc_particle_single_by_block_buf_t = psc_particle_buf<particle_single_by_block_t>;
#include "psc_particles_common.h"
#undef PTYPE

using mparticles_single_by_block_t = mparticles<psc_particle_single_by_block_buf_t, psc_mparticles_single_by_block>;

// can't do this as inline function since struct psc isn't known yet
#define particle_single_by_block_qni_div_mni(p) ({			\
      particle_single_by_block_real_t rv;				\
      rv = ppsc->kinds[p->kind].q / ppsc->kinds[p->kind].m;		\
      rv;								\
    })

#define particle_single_by_block_qni(p) ({				\
      particle_single_by_block_real_t rv;				\
      rv = ppsc->kinds[p->kind].q;					\
      rv;								\
    })

#define particle_single_by_block_mni(p) ({				\
      particle_single_by_block_real_t rv;				\
      rv = ppsc->kinds[p->kind].m;					\
      rv;								\
    })

#define particle_single_by_block_wni(p) ({				\
      particle_single_by_block_real_t rv;				\
      rv = p->qni_wni / ppsc->kinds[p->kind].q;				\
      rv;								\
    })

static inline particle_single_by_block_real_t
particle_single_by_block_qni_wni(particle_single_by_block_t *p)
{
  return p->qni_wni;
}

#endif
