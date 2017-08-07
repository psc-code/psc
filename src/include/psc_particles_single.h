
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc_particles_private.h"

#include "psc.h"

typedef float particle_single_real_t;

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

typedef struct psc_particle_single {
  particle_single_real_t xi, yi, zi;
  particle_single_real_t qni_wni;
  particle_single_real_t pxi, pyi, pzi;
  int kind;
} particle_single_t;

struct psc_particles_single {
  particle_single_t *particles;
  particle_single_t *particles_alt;
  int b_mx[3];
  int nr_blocks;
  particle_single_real_t b_dxi[3];
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int n_send;
  unsigned int n_part_save;
  bool need_reorder;
};

#define psc_particles_single(prts) mrc_to_subobj(prts, struct psc_particles_single)

static inline particle_single_t *
particles_single_get_one(struct psc_particles *prts, int n)
{
  assert(psc_particles_ops(prts) == &psc_particles_single_ops);
  return &psc_particles_single(prts)->particles[n];
}

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

#define particle_single_x(prt) ((prt)->xi)
#define particle_single_px(prt) ((prt)->pxi)

static inline int
particle_single_kind(particle_single_t *prt)
{
  return prt->kind;
}

static inline void
__calc_vxi(particle_single_real_t vxi[3], particle_single_t *part)
{
  particle_single_real_t root =
    1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
particle_single_get_relative_pos(particle_single_t *p, double xb[3],
				 particle_single_real_t xi[3])
{
  particle_single_real_t vxi[3];
  __calc_vxi(vxi, p);
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
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
particle_single_real_nint(particle_single_real_t x)
{
  return floorf(x + .5f); // FIXME use roundf()?
}

static inline int
particle_single_real_fint(particle_single_real_t x)
{
  return floorf(x);
}

static inline particle_single_real_t
particle_single_real_sqrt(particle_single_real_t x)
{
  return sqrtf(x);
}

static inline particle_single_real_t
particle_single_real_abs(particle_single_real_t x)
{
  return fabsf(x);
}

#endif
