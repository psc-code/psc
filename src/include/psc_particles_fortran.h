
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc_particles_private.h"

#define PTYPE PTYPE_FORTRAN
#include "psc_particle_buf_common.h"
#include "psc_particles_common.h"
#undef PTYPE

static inline int
particle_fortran_real_fint(particle_fortran_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
