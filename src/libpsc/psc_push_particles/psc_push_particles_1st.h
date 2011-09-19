
#ifndef PSC_PUSH_PARTICLES_1ST_H
#define PSC_PUSH_PARTICLES_1ST_H

#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

// switch between double and float in generic_c
// constants need to always be given like 1.5f
// they will propagated to double precision as necessary
// when actually doing double computations

#define CREAL 8

#if CREAL == 8

typedef double creal;
#define creal_abs(x) fabs(x)
#define creal_sqrt(x) sqrt(x)

#elif CREAL == 4

typedef float creal;
#define creal_abs(x) fabsf(x)
#define creal_sqrt(x) sqrtf(x)

#endif

void psc_push_particles_1st_push_xz(struct psc_push_particles *push,
				    mparticles_base_t *particles_base,
				    mfields_base_t *flds_base);

static inline int
nint(creal x)
{
  return (int)(x + 10.5f) - 10;
}

// like floor(), though returns int

static inline int
fint(creal x)
{
  return (int)(x + 10.f) - 10;
}

#endif
