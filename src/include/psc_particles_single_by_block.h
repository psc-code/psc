
#ifndef PSC_PARTICLE_SINGLE_BY_BLOCK_H
#define PSC_PARTICLE_SINGLE_BY_BLOCK_H

#include "psc_particles_private.h"

#include "psc.h"

typedef float particle_single_by_block_real_t;

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

typedef struct psc_particle_single_by_block {
  particle_single_by_block_real_t xi, yi, zi;
  particle_single_by_block_real_t qni_wni;
  particle_single_by_block_real_t pxi, pyi, pzi;
  int kind;
} particle_single_by_block_t;

struct psc_mparticles_single_by_block_patch {
  particle_single_by_block_t *prt_array;
  particle_single_by_block_t *prt_array_alt;
  int b_mx[3];
  int nr_blocks;
  particle_single_by_block_real_t b_dxi[3];
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int *b_off;
  bool need_reorder;
};

struct psc_mparticles_single_by_block {
  struct psc_mparticles_single_by_block_patch *patch;
};

#define psc_mparticles_single_by_block(prts) mrc_to_subobj(prts, struct psc_mparticles_single_by_block)

static inline particle_single_by_block_t *
particles_single_by_block_get_one(struct psc_particles *prts, int n)
{
  struct psc_mparticles *mprts = prts->mprts;
  int p = prts->p;
  
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_single_by_block_ops);
  return &psc_mparticles_single_by_block(mprts)->patch[p].prt_array[n];
}

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

#define particle_single_by_block_x(prt) ((prt)->xi)
#define particle_single_by_block_px(prt) ((prt)->pxi)

static inline particle_single_by_block_real_t
particle_single_by_block_qni_wni(particle_single_by_block_t *p)
{
  return p->qni_wni;
}

static inline int
particle_single_by_block_kind(particle_single_by_block_t *prt)
{
  return prt->kind;
}

static inline void
particle_single_by_block_calc_vxi(particle_single_by_block_real_t vxi[3],
				  particle_single_by_block_t *part)
{
  particle_single_by_block_real_t root =
    1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
particle_single_by_block_get_relative_pos(particle_single_by_block_t *p, double xb[3],
				 particle_single_by_block_real_t xi[3])
{
  particle_single_by_block_real_t vxi[3];
  particle_single_by_block_calc_vxi(vxi, p);
  particle_single_by_block_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
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
particle_single_by_block_real_nint(particle_single_by_block_real_t x)
{
  return floorf(x + .5f); // FIXME use roundf()?
}

static inline int
particle_single_by_block_real_fint(particle_single_by_block_real_t x)
{
  return floorf(x);
}

static inline particle_single_by_block_real_t
particle_single_by_block_real_sqrt(particle_single_by_block_real_t x)
{
  return sqrtf(x);
}

static inline particle_single_by_block_real_t
particle_single_by_block_real_abs(particle_single_by_block_real_t x)
{
  return fabsf(x);
}

#endif
