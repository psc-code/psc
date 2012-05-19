
#ifndef PSC_PARTICLE_DOUBLE_H
#define PSC_PARTICLE_DOUBLE_H

#include "psc_particles_private.h"

typedef float particle_double_real_t;

#define MPI_PARTICLES_DOUBLE_REAL MPI_FLOAT

typedef struct psc_particle_double {
  particle_double_real_t xi, yi, zi;
  particle_double_real_t qni_wni;
  particle_double_real_t pxi, pyi, pzi;
  int kind;
} particle_double_t;

struct psc_particles_double {
  particle_double_t *particles;
  int n_alloced;
};

#define psc_particles_double(prts) mrc_to_subobj(prts, struct psc_particles_double)

static inline particle_double_t *
particles_double_get_one(struct psc_particles *prts, int n)
{
  return &psc_particles_double(prts)->particles[n];
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

#define particle_double_wni(p) ({				\
      particle_double_real_t rv;				\
      rv = p->qni_wni / ppsc->kinds[p->kind].q;			\
      rv;							\
    })

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
particle_double_get_relative_pos(particle_double_t *p, double xb[3],
				 particle_double_real_t xi[3])
{
  xi[0] = p->xi;
  xi[1] = p->yi;
  xi[2] = p->zi;
}

static inline int
particle_double_real_nint(particle_double_real_t x)
{
  return (int)(x + 10.5f) - 10;
}

static inline int
particle_double_real_fint(particle_double_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
