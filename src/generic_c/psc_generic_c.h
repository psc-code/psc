
#include "psc.h"
#include "psc_fields_c.h"
#include "psc_particles_c.h"

// switch between double and float in generic_c
// constants need to always be given like 1.5f
// they will propagated to double precision as necessary
// when actually doing double computations

#define CREAL 8
#define FREAL 8

#if CREAL == FREAL
//#define USE_FF3
#endif

#if CREAL == 8

typedef double creal;
#define creal_abs(x) fabs(x)
#define creal_sqrt(x) sqrt(x)

#elif CREAL == 4

typedef float creal;
#define creal_abs(x) fabsf(x)
#define creal_sqrt(x) sqrtf(x)

#endif

// ----------------------------------------------------------------------
// generic C data structures

struct psc_genc {
  particle_c_t *part;
  int part_allocated;
  psc_fields_c_t pf;
};

void genc_push_part_xz();
void genc_push_part_yz_a();
void genc_push_part_yz_b();

#define F3(m, jx,jy,jz) F3_C(pf, m, jx,jy,jz)

static inline int
nint(creal x)
{
  return (int)(x + 10.5f) - 10;
}

