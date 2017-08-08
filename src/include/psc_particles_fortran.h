
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc_particles_private.h"

#define PTYPE PTYPE_FORTRAN
#include "psc_particles_common.h"
#undef PTYPE

static inline particle_fortran_t *
psc_mparticles_fortran_get_one(struct psc_mparticles *mprts, int p, int n)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_fortran_ops);
  return &psc_mparticles_fortran(mprts)->patch[p].prt_array[n];
}

static inline int
particle_fortran_real_fint(particle_fortran_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
