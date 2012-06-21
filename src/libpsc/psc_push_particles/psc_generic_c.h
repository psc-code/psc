
#ifndef PSC_GENERIC_C_H
#define PSC_GENERIC_C_H

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

void psc_push_particles_generic_c_push_a_y(struct psc_push_particles *push,
					   struct psc_particles *particles_base,
					   struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_z(struct psc_push_particles *push,
					   struct psc_particles *particles_base,
					   struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_xy(struct psc_push_particles *push,
					    struct psc_particles *particles_base,
					    struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_xz(struct psc_push_particles *push,
					    struct psc_particles *particles_base,
					    struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_yz(struct psc_push_particles *push,
					    struct psc_particles *particles_base,
					    struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_xyz(struct psc_push_particles *push,
					     struct psc_particles *particles_base,
					     struct psc_fields *flds_base);

void psc_push_particles_generic_c_calc_j_z(struct psc_push_particles *push,
					   struct psc_particles *prts_base,
					   struct psc_fields *flds_base);

void psc_push_particles_generic_c_push_yz_a(struct psc_push_particles *push,
					    mparticles_base_t *particles_base,
					    mfields_base_t *flds_base);
void psc_push_particles_generic_c_push_yz_b(struct psc_push_particles *push,
					    mparticles_base_t *particles_base,
					    mfields_base_t *flds_base);

static inline int
nint(creal x)
{
  return (int)(x + 10.5f) - 10;
}

#endif
