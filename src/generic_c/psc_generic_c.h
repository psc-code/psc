
#ifndef PSC_GENERIC_C_H
#define PSC_GENERIC_C_H

#include "psc.h"
#include "psc_fields_as_fortran.h"
#include "psc_particles_as_fortran.h"

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

void genc_push_part_xz();
void genc_push_part_yz();
void genc_push_part_z();
void genc_push_part_yz_a();
void genc_push_part_yz_b();

static inline int
nint(creal x)
{
  return (int)(x + 10.5f) - 10;
}

#endif
