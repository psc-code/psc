
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc_particles_private.h"

#include <psc.h>

#define PTYPE PTYPE_SINGLE
#include "psc_particle_buf_common.h"
#include "psc_particles_common.h"
#undef PTYPE

// can't do this as inline function since struct psc isn't known yet
#define particle_single_qni_div_mni(p) ({			\
      particle_single_real_t rv;				\
      rv = ppsc->kinds[p->kind].q / ppsc->kinds[p->kind].m;	\
      rv;							\
    })

#define particle_single_qni(p) ({				\
      particle_single_real_t rv;				\
      rv = ppsc->kinds[p->kind].q;				\
      rv;							\
    })

#define particle_single_mni(p) ({				\
      particle_single_real_t rv;				\
      rv = ppsc->kinds[p->kind].m;				\
      rv;							\
    })

#define particle_single_wni(p) ({				\
      particle_single_real_t rv;				\
      rv = p->qni_wni / ppsc->kinds[p->kind].q;			\
      rv;							\
    })

#define particle_single_qni_wni(prt) ((prt)->qni_wni)

static inline int
particle_single_kind(particle_single_t *prt)
{
  return prt->kind;
}

static inline int
particle_single_real_nint(particle_single_real_t x)
{
  return floorf(x + .5f); // FIXME use roundf()?
}

static inline particle_single_real_t
particle_single_real_sqrt(particle_single_real_t x)
{
  return sqrtf(x);
}

#endif
