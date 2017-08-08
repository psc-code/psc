
#ifndef PSC_PARTICLE_DOUBLE_H
#define PSC_PARTICLE_DOUBLE_H

#include "psc_particles_private.h"

#define PTYPE PTYPE_DOUBLE
#include "psc_particles_common.h"
#undef PTYPE

struct psc_mparticles_double {
  struct psc_mparticles_double_patch *patch;
};

#define psc_mparticles_double(prts) mrc_to_subobj(prts, struct psc_mparticles_double)

#include <math.h>
#include "psc.h"

static inline particle_double_t *
psc_mparticles_double_get_one(struct psc_mparticles *mprts, int p, int n)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_double_ops);
  return &psc_mparticles_double(mprts)->patch[p].prt_array[n];
}

// can't do this as inline function since struct psc isn't known yet
#define particle_double_qni_div_mni(p) ({			\
      particle_double_real_t rv;				\
      rv = ppsc->kinds[p->kind].q / ppsc->kinds[p->kind].m;	\
      rv;							\
    })

#define particle_double_qni(p) ({				\
      particle_double_real_t rv;				\
      rv = ppsc->kinds[p->kind].q;				\
      rv;							\
    })

#define particle_double_mni(p) ({				\
      particle_double_real_t rv;				\
      rv = ppsc->kinds[p->kind].m;				\
      rv;							\
    })

#define particle_double_wni(p) ({				\
      particle_double_real_t rv;				\
      rv = p->qni_wni / ppsc->kinds[p->kind].q;			\
      rv;							\
    })

#define particle_double_x(prt) ((prt)->xi)
#define particle_double_px(prt) ((prt)->pxi)

static inline particle_double_real_t
particle_double_qni_wni(particle_double_t *p)
{
  return p->qni_wni;
}

static inline int
particle_double_kind(particle_double_t *prt)
{
  return prt->kind;
}

static inline void
particle_double_calc_vxi(particle_double_real_t vxi[3], particle_double_t *part)
{
  particle_double_real_t root =
    1.f / sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
particle_double_get_relative_pos(particle_double_t *p, double xb[3],
				 particle_double_real_t xi[3])
{
  particle_double_real_t vxi[3];
  particle_double_calc_vxi(vxi, p);
  particle_double_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  xi[0] = p->xi - dth[0] * vxi[0];
  xi[1] = p->yi - dth[1] * vxi[1];
  xi[2] = p->zi - dth[2] * vxi[2];
}

static inline int
particle_double_real_nint(particle_double_real_t x)
{
  return floor(x + .5);
}

static inline int
particle_double_real_fint(particle_double_real_t x)
{
  return floor(x);
}

static inline particle_double_real_t
particle_double_real_sqrt(particle_double_real_t x)
{
  return sqrt(x);
}

static inline particle_double_real_t
particle_double_real_abs(particle_double_real_t x)
{
  return fabs(x);
}

#endif
